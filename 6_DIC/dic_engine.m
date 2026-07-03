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
%
%   Output struct R:
%     R.U, R.V   displacements, px (x right / y down), NaN where invalid
%     R.Q        ZNCC quality per point
%     R.NIT      IC-GN iteration count per point
%     R.gx, R.gy grid coordinates (columns / rows of the image)
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
    nY = numel(gy); nX = numel(gx);

    % mask of "black" subsets (the hole): max brightness in the window
    refmax = movmax(movmax(ref, P.subset_size, 1), P.subset_size, 2);
    speckle = refmax(gy, gx) >= P.black_thr;

    % spline interpolant of the deformed frame (built once)
    Fi = griddedInterpolant({1:Himg, 1:Wimg}, cur, 'spline', 'nearest');

    % subpixel basis: local subset coordinates
    [xi, eta] = meshgrid(-half:half, -half:half);  % xi = X, eta = Y
    xi = xi(:); eta = eta(:);

    U = nan(nY, nX); V = nan(nY, nX); Q = nan(nY, nX); NIT = nan(nY, nX);
    Pw = nan(nY, nX, 6);   % parameters [u ux uy v vx vy] per point

    % --- seed: the speckle node closest to the image center ---
    [Jg, Ig] = meshgrid(1:nX, 1:nY);
    d2 = (Ig - (nY+1)/2).^2 + (Jg - (nX+1)/2).^2;
    d2(~speckle) = inf;
    [~, kseed] = min(d2(:));
    [si, sj] = ind2sub([nY, nX], kseed);

    [du0, dv0] = int_search(ref, cur, gy(si), gx(sj), half, P.search_range);
    [p, q, nit, ok] = icgn_point(ref, Fi, gy(si), gx(sj), ...
        [du0 0 0 dv0 0 0], half, xi, eta, P);
    if ~ok, error('DIC: seed point did not converge - increase search_range.'); end
    Pw(si, sj, :) = p;
    U(si, sj) = p(1); V(si, sj) = p(4); Q(si, sj) = q; NIT(si, sj) = nit;

    % --- RG-DIC: front grows from the best ZNCC (binary heap) ---
    % A node is finalized ONLY on success: a node that failed from one
    % neighbor can still converge from another (no dead holes), and the
    % heap keeps dense grids at O(N log N).
    cap = nY * nX;
    hq = zeros(cap, 1); hi = zeros(cap, 1); hj = zeros(cap, 1); hn = 0;
    [hq, hi, hj, hn] = heap_push(hq, hi, hj, hn, q, si, sj);
    done = false(nY, nX);
    done(si, sj) = true;
    nbr = [1 0; -1 0; 0 1; 0 -1];
    while hn > 0
        [ci, cj, hq, hi, hj, hn] = heap_pop(hq, hi, hj, hn);
        p0 = squeeze(Pw(ci, cj, :)).';
        for n = 1:4
            ii = ci + nbr(n,1); jj = cj + nbr(n,2);
            if ii < 1 || ii > nY || jj < 1 || jj > nX, continue; end
            if done(ii, jj) || ~speckle(ii, jj), continue; end
            [p, q, nit, ok] = icgn_point(ref, Fi, gy(ii), gx(jj), transfer(p0, nbr(n,:), P.grid_step), half, xi, eta, P);
            if ~ok || q < P.zncc_min, continue; end
            done(ii, jj) = true;
            Pw(ii, jj, :) = p;
            U(ii, jj) = p(1); V(ii, jj) = p(4); Q(ii, jj) = q; NIT(ii, jj) = nit;
            [hq, hi, hj, hn] = heap_push(hq, hi, hj, hn, q, ii, jj);
        end
    end
    R = struct('U', U, 'V', V, 'Q', Q, 'NIT', NIT, 'gx', gx, 'gy', gy);
end

function [du, dv] = int_search(ref, cur, y0, x0, half, rng_px)
% Integer ZNCC search around (x0, y0) - used only for the seed point.
    f = ref(y0-half:y0+half, x0-half:x0+half);
    fz = f(:) - mean(f(:)); fn = norm(fz);
    [Himg, Wimg] = size(cur);
    best = -2; du = 0; dv = 0;
    for dy = -rng_px:rng_px
        for dx = -rng_px:rng_px
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

function [p, zncc, nit, ok] = icgn_point(ref, Fi, y0, x0, p0, half, xi, eta, P)
% IC-GN at one point: ZNSSD criterion, first-order shape function
% p = [u ux uy v vx vy]; inverse-compositional warp update
% W <- W * inv(W(dp)) (Baker-Matthews), residual e = I - T.
    p = p0; zncc = -2; nit = 0; ok = false;
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
        if c < P.conv_tol, break; end
    end
    p = [W(1,3), W(1,1)-1, W(1,2), W(2,3), W(2,1), W(2,2)-1];
    zncc = (fz' * gz) / (fn * gn);
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
