function R = dic_engine(ref, cur, P)
% DIC_ENGINE  Subpixel 2D DIC correlation engine: ZNSSD + IC-GN + RG-DIC.
%
%   R = DIC_ENGINE(ref, cur, P) registers the deformed frame `cur` against
%   the reference frame `ref` (both grayscale, double, in [0..1]) and
%   returns the displacement fields on a regular grid. The engine knows
%   nothing about where the frames come from or how the fields will be
%   used; any script (dic_icgn.m, the Challenge 2.0 benchmark, ...) can
%   call it on arbitrary inputs.
%
%   Parameters (struct P):
%     P.subset_size   correlation subset, px (odd)
%     P.grid_step     grid step, px
%     P.search_range  integer search range, px (seed point only)
%     P.black_thr     "black region" threshold (masked-out areas)
%     P.zncc_min      minimum quality; points below are rejected
%     P.max_iter      IC-GN iteration cap
%     P.conv_tol      convergence tolerance on ||dp||, px
%     P.roi           optional [x y w h] region of interest, px; the grid
%                     is restricted to nodes inside it (default: full frame)
%     P.seed_shift    optional [sx sy] coarse prior shift, px: the integer
%                     seed search of +/-search_range is CENTERED on it, so
%                     large rigid-body motion (hundreds of px) can be
%                     handled by an external coarse estimator without
%                     enlarging the search window (default: [0 0])
%     P.seed_xy       optional [x y] preferred seed location, px: the seed
%                     becomes the speckle grid node nearest to that point
%                     instead of the grid center - callers can place it on
%                     the best-textured spot (default: grid center)
%     P.mask          optional full-frame boolean the SAME size as `ref`:
%                     grid nodes whose pixel is false are excluded (kept
%                     NaN), on top of P.roi and the black-region mask. Lets
%                     an arbitrary AOI - e.g. a polygon that cuts out
%                     fixtures/grips/glare - be passed in without changing
%                     the rectangular P.roi. Absent/empty = no extra mask
%                     (default), so existing callers are unaffected.
%     P.multi_seed    restart the propagation front in unreached speckle
%                     regions (default true). When a correlation barrier
%                     crosses the whole grid (weld bead, tape, glare band),
%                     a single front dies there and everything beyond stays
%                     NaN even though it is perfectly trackable; re-seeding
%                     recovers those regions. On barrier-free data no
%                     re-seed ever triggers and the result is bit-identical
%                     to the single-seed engine. false = classical RG-DIC.
%     P.max_seeds     cap on the number of successful seeds, incl. the
%                     primary one. Default: as many as the grid could
%                     legitimately need, i.e. how many regions of the
%                     minimum seedable size fit in the speckle area. A
%                     seed is only ever attempted while that much of the
%                     grid is still unreached, so this is only a runaway
%                     guard (see the loop).
%     P.max_seed_fails  seed attempts that may find nothing before the
%                     engine stops looking (default 200). Separate from
%                     max_seeds: a successful seed fills a region, a
%                     failed one is only cost, so only failures are
%                     bounded.
%
%   Output struct R:
%     R.U, R.V   displacements, px (x right / y down), NaN where invalid
%     R.Q        ZNCC quality per point
%     R.NIT      IC-GN iteration count per point
%     R.CONV     true where IC-GN reached P.conv_tol; false where the point
%                was kept after exhausting P.max_iter (its Q still passed
%                zncc_min, but the solve itself did not converge)
%     R.nonconv  fraction of accepted points with R.CONV false
%     R.SU, R.SV per-point standard uncertainty of U and V, px (1 sigma):
%                the Gauss-Newton covariance at the solution,
%                Cov(p) = s^2 * inv(H) with s^2 = |e|^2 / (N - 6) the
%                residual variance of the subset (Wang, Pan et al. 2009;
%                the same estimate DICe and VIC-3D report as "sigma").
%                It is the RANDOM part of the error - image noise seen
%                through the local texture - not the interpolation bias
%                (that one is measured separately by bias_scurve.m).
%     R.gx, R.gy grid coordinates (columns / rows of the image)
%     R.nseeds   number of seeds actually used (1 = single front sufficed)
%     R.nfail    seed attempts that found nothing
%     R.seed_capped  true when max_seeds ran out rather than the grid: the
%                    NaN left behind is a budget, not untrackable data, and
%                    nothing else in the output distinguishes the two
%
%   Method: ZNSSD criterion (invariant to linear illumination changes);
%   inverse-compositional Gauss-Newton (Baker & Matthews) with a
%   first-order (affine) subset shape function and a Hessian assembled
%   once per point; spline interpolation of the deformed frame;
%   reliability-guided propagation (RG-DIC, Pan) from a single seed with
%   a binary-heap front and first-order initial-guess transfer
%   (Pan, Li & Tong 2013). No toolboxes and no GPU required.
    half = floor(P.subset_size / 2);
    [Himg, Wimg] = size(ref);
    margin = half + P.search_range + 2;
    gy = (margin+1):P.grid_step:(Himg-margin);   % rows (Y)
    gx = (margin+1):P.grid_step:(Wimg-margin);   % columns (X)
    if isfield(P, 'roi') && ~isempty(P.roi)      % optional region [x y w h]
        gy = gy(gy >= P.roi(2) & gy <= P.roi(2) + P.roi(4));
        gx = gx(gx >= P.roi(1) & gx <= P.roi(1) + P.roi(3));
        if isempty(gy) || isempty(gx)
            error(['DIC: the ROI leaves no grid nodes - enlarge the ROI ' ...
                   'or reduce subset/step.']);
        end
    end
    nY = numel(gy); nX = numel(gx);

    % mask of "black" subsets (the hole): max brightness in the window
    refmax = movmax(movmax(ref, P.subset_size, 1), P.subset_size, 2);
    speckle = refmax(gy, gx) >= P.black_thr;

    % optional user AOI mask: a full-frame boolean sampled at the grid nodes
    % and AND-ed into the speckle mask. Nodes it excludes are never seeded or
    % grown, so they stay NaN - exactly like the black-region mask. This is
    % the engine half of the polygon-AOI feature (the interactive editor is
    % UI); any caller can restrict the grid to an arbitrary region without
    % touching the rectangular P.roi.
    if isfield(P, 'mask') && ~isempty(P.mask)
        if ~isequal(size(P.mask), [Himg, Wimg])
            error(['DIC: P.mask must be the same size as the reference ' ...
                   'frame (%dx%d); got %dx%d.'], Himg, Wimg, ...
                   size(P.mask, 1), size(P.mask, 2));
        end
        speckle = speckle & logical(P.mask(gy, gx));
        if ~any(speckle(:))
            error(['DIC: P.mask (with P.roi) leaves no valid grid nodes - ' ...
                   'widen the mask/ROI or reduce subset/step.']);
        end
    end

    % spline interpolant of the deformed frame (built once)
    Fi = griddedInterpolant({1:Himg, 1:Wimg}, cur, 'spline', 'nearest');

    % subpixel basis: local subset coordinates
    [xi, eta] = meshgrid(-half:half, -half:half);  % xi = X, eta = Y
    xi = xi(:); eta = eta(:);

    U = nan(nY, nX); V = nan(nY, nX); Q = nan(nY, nX); NIT = nan(nY, nX);
    CONV = false(nY, nX);
    SU = nan(nY, nX); SV = nan(nY, nX);   % 1-sigma uncertainty per point, px
    Pw = nan(nY, nX, 6);   % parameters [u ux uy v vx vy] per point

    % --- seed: the speckle node closest to the requested point (default:
    %     the grid center) ---
    [Jg, Ig] = meshgrid(1:nX, 1:nY);
    if isfield(P, 'seed_xy') && ~isempty(P.seed_xy)
        d2 = (gy(Ig) - P.seed_xy(2)).^2 + (gx(Jg) - P.seed_xy(1)).^2;
    else
        d2 = (Ig - (nY+1)/2).^2 + (Jg - (nX+1)/2).^2;
    end
    d2(~speckle) = inf;
    [~, kseed] = min(d2(:));
    [si, sj] = ind2sub([nY, nX], kseed);

    if isfield(P, 'seed_shift') && ~isempty(P.seed_shift)
        s0 = round(double(P.seed_shift(:)'));    % [sx sy] coarse prior
    else
        s0 = [0 0];
    end

    % --- RG-DIC with multi-seed restart (binary heap front) ---
    % A node is finalized ONLY on success: a node that failed from one
    % neighbor can still converge from another (no dead holes), and the
    % heap keeps dense grids at O(N log N).
    % One seed + one front is classical RG-DIC. When a correlation barrier
    % crosses the WHOLE grid (weld bead, tape, glare band), the front dies
    % at it and every node beyond stays NaN even though it is perfectly
    % trackable. So after the front exhausts, the engine re-seeds at the
    % best-textured node of the remaining unreached speckle area (same
    % integer search around the same coarse prior) and grows a new front.
    % Re-seeds only ever fill nodes that are still NaN, so on barrier-free
    % data (no unreached nodes) the result is bit-identical to the
    % single-seed engine. P.multi_seed = false restores the old behavior.
    cap = nY * nX;
    hq = zeros(cap, 1); hi = zeros(cap, 1); hj = zeros(cap, 1); hn = 0;
    done  = false(nY, nX);
    tried = false(nY, nX);   % nodes rejected as re-seed candidates
    nbr = [1 0; -1 0; 0 1; 0 -1];
    multi = ~isfield(P, 'multi_seed') || isempty(P.multi_seed) || P.multi_seed;
    minregion = max(10, round(0.005 * nnz(speckle)));   % ignore tiny islands
    % When to stop re-seeding. A fixed cap (was 8) does not work: number of
    % fronts needed depends on grid sampling. DIC 3D dataset, a specimen speckled in two
    % zones, cross-camera match at grid step 8 / subset 15:
    %
    %     seeds allowed      8      59     201    1030 (unbounded)
    %     far zone         DEAD    full    full    full
    %     near zone, valid nodes per row
    %                       72      72      72      74
    %     seconds          7.8    13.5    18.9    47.1
    %
    % At 8 the far zone (7800 nodes, a fifth of the grid) was never seeded
    % from the third state on. At the default step 20 a few seeds are
    % enough, so it did not show up there.
    %
    % Past ~60 seeds the tail is scattered fragments, but a "stop when
    % yield drops" rule does not work here. Seeds are picked by TEXTURE, not by region size, so a run of small
    % fragments can precede a large unreached zone: with a stop after five
    % consecutive low-yield seeds, state 21-0547 of this same series halted
    % at nine seeds having never reached the far zone, and the comparison
    % fell straight back to 39 % of the specimen. A yield rule needs the
    % seeds ordered by region size to be sound, which costs a connected
    % component pass per seed.
    %
    % So: seed while a region worth seeding remains, up to as many regions
    % of the minimum size as fit in the grid; cost is bounded by the failure
    % budget below. R.seed_capped reports if the bound was reached.
    maxseeds = max(8, ceil(nnz(speckle) / minregion));
    if isfield(P, 'max_seeds') && ~isempty(P.max_seeds), maxseeds = P.max_seeds; end
    % Failure budget is separate: a failed seed is a full integer search
    % over +/-search_range for nothing. Tied to maxseeds, raising the seed
    % budget would also raise the worst case.
    MAXFAIL = 200;
    if isfield(P, 'max_seed_fails') && ~isempty(P.max_seed_fails)
        MAXFAIL = P.max_seed_fails;
    end
    nseeds = 0; nfail = 0;
    while true
        % ---- pick a seed node ----
        if nseeds > 0            % primary (si,sj) already chosen above
            if ~multi || nseeds >= maxseeds, break; end
            cand = speckle & ~done & ~tried;
            if nnz(cand) < minregion, break; end
            [si, sj] = best_texture_node(ref, gy, gx, cand, half);
        end
        % ---- solve the seed ----
        [du0, dv0] = int_search(ref, cur, gy(si), gx(sj), half, P.search_range, s0);
        [p, q, nit, ok, cv, sg] = icgn_point(ref, Fi, gy(si), gx(sj), ...
            [du0 0 0 dv0 0 0], half, xi, eta, P);
        if nseeds == 0 && ~ok
            error('DIC: seed point did not converge - increase search_range.');
        end
        if ~ok || (nseeds > 0 && q < P.zncc_min)
            % bad seed spot (e.g. ON the barrier): block its 3x3 ring too
            tried(max(1,si-1):min(nY,si+1), max(1,sj-1):min(nX,sj+1)) = true;
            nfail = nfail + 1;
            if nfail > MAXFAIL, break; end       % bad-region bailout
            continue;
        end
        nseeds = nseeds + 1;
        done(si, sj) = true;
        Pw(si, sj, :) = p;
        U(si, sj) = p(1); V(si, sj) = p(4); Q(si, sj) = q; NIT(si, sj) = nit;
        CONV(si, sj) = cv; SU(si, sj) = sg(1); SV(si, sj) = sg(2);
        [hq, hi, hj, hn] = heap_push(hq, hi, hj, hn, q, si, sj);
        % ---- grow the front (classical RG-DIC, unchanged) ----
        while hn > 0
            [ci, cj, hq, hi, hj, hn] = heap_pop(hq, hi, hj, hn);
            p0 = squeeze(Pw(ci, cj, :)).';
            for n = 1:4
                ii = ci + nbr(n,1); jj = cj + nbr(n,2);
                if ii < 1 || ii > nY || jj < 1 || jj > nX, continue; end
                if done(ii, jj) || ~speckle(ii, jj), continue; end
                [p, q, nit, ok, cv, sg] = icgn_point(ref, Fi, gy(ii), gx(jj), transfer(p0, nbr(n,:), P.grid_step), half, xi, eta, P);
                if ~ok || q < P.zncc_min, continue; end
                done(ii, jj) = true;
                Pw(ii, jj, :) = p;
                U(ii, jj) = p(1); V(ii, jj) = p(4); Q(ii, jj) = q; NIT(ii, jj) = nit;
                CONV(ii, jj) = cv; SU(ii, jj) = sg(1); SV(ii, jj) = sg(2);
                [hq, hi, hj, hn] = heap_push(hq, hi, hj, hn, q, ii, jj);
            end
        end
    end
    % Speckle the propagation never got to. Distinct from "failed to
    % correlate": these nodes were never attempted at all.
    %
    % Two numbers, the total alone is not enough. DIC 3D dataset, one
    % cross-camera pair: total unreached 46.6 % with a zone lost, 13.8 % at
    % the seed cap, 9.2 % unbounded - no threshold on it works. The largest
    % connected missing region separates them: a lost zone is one big area,
    % a healthy run leaves small scattered fragments.
    nsp = max(1, nnz(speckle));
    miss = speckle & ~done;
    unreached = nnz(miss) / nsp;
    unreached_max = largest_region(miss) / nsp;
    % Accepted points whose IC-GN solve hit max_iter instead of conv_tol.
    % They pass the zncc_min gate on SIMILARITY, not on convergence of the
    % SOLUTION, so they are kept in the field; the mask and its fraction
    % are returned next to NIT.
    nonconv = nnz(done & ~CONV) / max(1, nnz(done));
    R = struct('U', U, 'V', V, 'Q', Q, 'NIT', NIT, 'CONV', CONV, ...
        'SU', SU, 'SV', SV, 'gx', gx, 'gy', gy, ...
        'nseeds', nseeds, 'nfail', nfail, 'seed_capped', nseeds >= maxseeds, ...
        'unreached', unreached, 'unreached_max', unreached_max, ...
        'nonconv', nonconv);
end

function n = largest_region(M)
% Size of the largest 4-connected true region of a logical grid. An
% explicit stack flood fill: BWLABEL and BWCONNCOMP are toolbox functions
% and this engine is base MATLAB throughout.
    n = 0;
    if ~any(M(:)), return; end
    [nY, nX] = size(M);
    seen = false(nY, nX);
    stack = zeros(numel(M), 2);
    for a = 1:nY
        for b = 1:nX
            if ~M(a, b) || seen(a, b), continue; end
            sp = 1; stack(1, :) = [a b]; seen(a, b) = true; sz = 0;
            while sp > 0
                i = stack(sp, 1); j = stack(sp, 2); sp = sp - 1;
                sz = sz + 1;
                for d = [1 0; -1 0; 0 1; 0 -1].'
                    ii = i + d(1); jj = j + d(2);
                    if ii < 1 || ii > nY || jj < 1 || jj > nX, continue; end
                    if M(ii, jj) && ~seen(ii, jj)
                        seen(ii, jj) = true;
                        sp = sp + 1; stack(sp, :) = [ii jj];
                    end
                end
            end
            if sz > n, n = sz; end
        end
    end
end

function [si, sj] = best_texture_node(ref, gy, gx, cand, half)
% Best re-seed candidate among the unreached speckle nodes: the node whose
% reference subset has the highest intensity std (strongest texture) - the
% safest place to restart the propagation front. A plain loop over the
% candidate nodes; each visit costs one subset std, so even thousands of
% candidates are cheap next to a single IC-GN solve.
    idx = find(cand);
    best = -inf; si = 0; sj = 0;
    nY = size(cand, 1);
    for t = 1:numel(idx)
        i = mod(idx(t) - 1, nY) + 1;
        j = floor((idx(t) - 1) / nY) + 1;
        f = ref(gy(i)-half:gy(i)+half, gx(j)-half:gx(j)+half);
        v = std(f(:));
        if v > best, best = v; si = i; sj = j; end
    end
end

function [du, dv] = int_search(ref, cur, y0, x0, half, rng_px, s0)
% Integer ZNCC search around (x0 + s0(1), y0 + s0(2)) - seed point only.
% s0 is an optional coarse prior shift [sx sy]; the +/-rng_px window is
% centered on it (default [0 0], the historical behavior).
    if nargin < 7, s0 = [0 0]; end
    f = ref(y0-half:y0+half, x0-half:x0+half);
    fz = f(:) - mean(f(:)); fn = norm(fz);
    [Himg, Wimg] = size(cur);
    best = -2; du = s0(1); dv = s0(2);
    for dy = s0(2)-rng_px : s0(2)+rng_px
        for dx = s0(1)-rng_px : s0(1)+rng_px
            ys = y0 + dy; xs = x0 + dx;
            if ys-half < 1 || ys+half > Himg || xs-half < 1 || xs+half > Wimg
                continue;
            end
            g = cur(ys-half:ys+half, xs-half:xs+half);
            gz = g(:) - mean(g(:)); gn = norm(gz);
            if fn < eps || gn < eps, continue; end
            c = (fz' * gz) / (fn * gn);
            if c > best, best = c; du = dx; dv = dy; end
        end
    end
end

function [p, zncc, nit, ok, conv, sig] = icgn_point(ref, Fi, y0, x0, p0, half, xi, eta, P)
% IC-GN at one point: ZNSSD criterion, first-order shape function
% p = [u ux uy v vx vy]; inverse-compositional warp update
% W <- W * inv(W(dp)) (Baker-Matthews), residual e = I - T.
% ok   = a solution was produced (subset stayed in frame, Hessian sane);
% conv = the update norm fell below P.conv_tol before P.max_iter ran out.
% ok without conv is a point the caller may keep but must not call converged.
% sig  = [sigma_u sigma_v], px: standard uncertainty of the displacement
%        from the Gauss-Newton covariance s^2 * inv(H), where s^2 is the
%        residual variance of the last iteration (the residual of the
%        solution, since the last update is below conv_tol). NaN when
%        no solution was produced.
    p = p0; zncc = -2; nit = 0; ok = false; conv = false; sig = [NaN NaN];
    f = ref(y0-half:y0+half, x0-half:x0+half);
    fz = f(:) - mean(f(:)); fn = norm(fz);
    if fn < eps, return; end

    % reference subset gradients (5-point central difference, O(h^4)) - computed once
    wp = ref(y0-half-2:y0+half+2, x0-half-2:x0+half+2);     % 2-px pad for the stencil
    gxr = (wp(3:end-2,1:end-4) - 8*wp(3:end-2,2:end-3) + 8*wp(3:end-2,4:end-1) - wp(3:end-2,5:end)) / 12;   % d/dx, [1 -8 0 8 -1]/12
    gyr = (wp(1:end-4,3:end-2) - 8*wp(2:end-3,3:end-2) + 8*wp(4:end-1,3:end-2) - wp(5:end,3:end-2)) / 12;   % d/dy, [1 -8 0 8 -1]/12

    J = [gxr(:), gxr(:).*xi, gxr(:).*eta, gyr(:), gyr(:).*xi, gyr(:).*eta];
    H = J' * J;
    if rcond(H) < 1e-12, return; end

    W = warp_mat(p);
    [Himg, Wimg] = size(ref);
    gz = []; gn = 0;
    for nit = 1:P.max_iter
        xw = x0 + W(1,1)*xi + W(1,2)*eta + W(1,3);
        yw = y0 + W(2,1)*xi + W(2,2)*eta + W(2,3);
        if min(xw) < 2 || max(xw) > Wimg-1 || min(yw) < 2 || max(yw) > Himg-1
            return;   % subset left the frame
        end
        g = Fi(yw, xw);
        gz = g - mean(g); gn = norm(gz);
        if gn < eps, return; end
        e  = (fn/gn) * gz - fz;                  % I - T (ZNSSD normalization)
        dp = H \ (J' * e);
        W  = W / warp_mat(dp);                   % inverse-compositional step
        c  = sqrt(dp(1)^2 + (dp(2)*half)^2 + (dp(3)*half)^2 + ...
                  dp(4)^2 + (dp(5)*half)^2 + (dp(6)*half)^2);
        if c < P.conv_tol, conv = true; break; end
    end
    p = [W(1,3), W(1,1)-1, W(1,2), W(2,3), W(2,1), W(2,2)-1];
    zncc = (fz' * gz) / (fn * gn);
    % Per-point uncertainty: the residual of both frames' noise seen through
    % the subset texture. 6 = parameters of the shape function.
    s2 = (e' * e) / max(1, numel(e) - 6);
    % Noise in the reference gradients inflates H (Pan 2009): the 5-point
    % stencil passes 130/144 of the pixel-noise variance into each gradient,
    % and both frames carry about half of s2 each. Remove that expected
    % contribution so sigma is not underestimated on noisy images.
    vn = 0.5 * s2 * 130 / 144;
    B = [ones(numel(xi),1), xi(:), eta(:)];
    D = vn * (B' * B);
    Hc = H - blkdiag(D, D);
    % The correction must leave Hc positive definite. rcond and a positive
    % diagonal do not guarantee that: with an indefinite Hc, Hc \ eye(6)
    % goes to LU, Hi(1,1) < 0 and sig becomes COMPLEX (seen on points with
    % ZNCC ~0.85).
    [~, notpd] = chol(Hc);
    if notpd || rcond(Hc) < 1e-12 || any(diag(Hc) <= 0), Hc = H; end
    Hi = Hc \ eye(6);
    sig = sqrt(s2 * [Hi(1,1), Hi(4,4)]);
    ok = true;
end

function W = warp_mat(p)
% 3x3 affine warp matrix from parameters [u ux uy v vx vy].
    W = [1+p(2), p(3), p(1);  p(5), 1+p(6), p(4);  0, 0, 1];
end

function p0 = transfer(p, d, step)
% First-order (Taylor) initial-guess transfer to a 4-neighbor (Pan, Li &
% Tong 2013, Eq. 8, generalized to the signed neighbor direction): the four
% displacement gradients carry over unchanged, and the neighbor's
% translation is predicted one grid step ahead from those gradients. This
% gives the IC-GN solve a much closer warm start than copying p verbatim,
% which matters most under large strain/rotation. p = [u ux uy v vx vy];
% d = [di dj] neighbor offset in grid NODES (di rows/Y, dj cols/X);
% step = grid step in px.
    p0 = p;                              % gradients (p0(2,3,5,6)) unchanged
    dx = d(2) * step;                    % offset in x (columns)
    dy = d(1) * step;                    % offset in y (rows)
    p0(1) = p(1) + p(2)*dx + p(3)*dy;    % u <- u + ux*dx + uy*dy
    p0(4) = p(4) + p(5)*dx + p(6)*dy;    % v <- v + vx*dx + vy*dy
end

function [hq, hi, hj, hn] = heap_push(hq, hi, hj, hn, q, i, j)
% Insert into a binary max-heap (keyed by q). Arrays are preallocated.
    hn = hn + 1;
    hq(hn) = q; hi(hn) = i; hj(hn) = j;
    k = hn;
    while k > 1
        pk = floor(k/2);
        if hq(pk) >= hq(k), break; end
        [hq(pk), hq(k)] = deal(hq(k), hq(pk));
        [hi(pk), hi(k)] = deal(hi(k), hi(pk));
        [hj(pk), hj(k)] = deal(hj(k), hj(pk));
        k = pk;
    end
end

function [i, j, hq, hi, hj, hn] = heap_pop(hq, hi, hj, hn)
% Extract the maximum from the binary heap.
    i = hi(1); j = hj(1);
    hq(1) = hq(hn); hi(1) = hi(hn); hj(1) = hj(hn);
    hn = hn - 1;
    k = 1;
    while true
        l = 2*k; r = l + 1; m = k;
        if l <= hn && hq(l) > hq(m), m = l; end
        if r <= hn && hq(r) > hq(m), m = r; end
        if m == k, break; end
        [hq(m), hq(k)] = deal(hq(k), hq(m));
        [hi(m), hi(k)] = deal(hi(k), hi(m));
        [hj(m), hj(k)] = deal(hj(k), hj(m));
        k = m;
    end
end
