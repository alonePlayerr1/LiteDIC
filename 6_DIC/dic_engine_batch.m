function R = dic_engine_batch(ref, cur, P)
% DIC_ENGINE_BATCH  Batched IC-GN for dense sweeps on smooth, small fields.
%
%   R = DIC_ENGINE_BATCH(ref, cur, P) computes the same fields as
%   DIC_ENGINE (same grid, same ZNSSD + IC-GN mathematics, same output
%   struct), but solves ALL grid points simultaneously as one vectorized
%   batch instead of point-by-point with reliability-guided propagation.
%   Every point is initialized from a single integer search at the image
%   center, so the method is valid when the displacement field stays
%   within the IC-GN convergence basin of that common initial guess
%   (a few px) - true for benchmark sets like the DIC Challenge Stars and
%   for self-test sweeps. For general fields (large or discontinuous
%   motion) use DIC_ENGINE, whose RG-DIC propagation warm-starts each
%   point from a converged neighbor.
%
%   Why it exists: the per-point engine spends most of its wall time in
%   interpreter overhead (one function call and one small interpolation
%   per point per iteration). Batching evaluates the spline interpolant
%   for millions of pixels per call and updates all warps with paged
%   matrix operations, so MATLAB's built-in multithreading uses every
%   CPU core - no toolboxes required. Speedups of two orders of
%   magnitude are typical on dense (grid_step = 1) sweeps.
%
%   Numerical equivalence with DIC_ENGINE (identical criterion, gradient
%   stencil, interpolant, update rule and convergence test; both converge
%   to the same optimum within P.conv_tol) is asserted by the callers:
%   RUN_TESTS (synthetic shift) and challenge2_benchmark.m (real images).
%
%   Inputs/outputs: identical to DIC_ENGINE.

    half = floor(P.subset_size / 2);
    [Himg, Wimg] = size(ref);
    margin = half + P.search_range + 2;
    gy = (margin+1):P.grid_step:(Himg-margin);   % rows (Y)
    gx = (margin+1):P.grid_step:(Wimg-margin);   % columns (X)
    nY = numel(gy); nX = numel(gx);

    % mask of "black" subsets (identical to dic_engine)
    refmax = movmax(movmax(ref, P.subset_size, 1), P.subset_size, 2);
    speckle = refmax(gy, gx) >= P.black_thr;

    % spline interpolant of the deformed frame (built once)
    Fi = griddedInterpolant({1:Himg, 1:Wimg}, cur, 'spline', 'nearest');

    % full-image reference gradients, 5-point stencil [1 -8 0 8 -1]/12
    % (correlation form; equals the per-subset slicing in dic_engine
    % because every subset pixel is >= 2 px away from the image border)
    kx = [-1 8 0 -8 1] / 12;                     % conv kernel = flipped stencil
    gxf = conv2(ref, kx,  'same');               % d/dx
    gyf = conv2(ref, kx.', 'same');              % d/dy

    % subset-local coordinates
    [xi, eta] = meshgrid(-half:half, -half:half);
    xi = xi(:); eta = eta(:);                    % S x 1
    S = numel(xi);

    % single integer-search seed at the speckle node closest to the center;
    % its offset initializes EVERY point (see the validity note above)
    [Jg, Ig] = meshgrid(1:nX, 1:nY);
    d2 = (Ig - (nY+1)/2).^2 + (Jg - (nX+1)/2).^2;
    d2(~speckle) = inf;
    [~, kseed] = min(d2(:));
    [si, sj] = ind2sub([nY, nX], kseed);
    [du0, dv0] = int_search(ref, cur, gy(si), gx(sj), half, P.search_range);

    % point list (speckle nodes only), image coordinates of the nodes
    [Ilist, Jlist] = find(speckle);
    y0 = gy(Ilist); x0 = gx(Jlist);              % 1 x M
    M = numel(y0);

    U = nan(nY, nX); V = nan(nY, nX); Q = nan(nY, nX); NIT = nan(nY, nX);

    % ---- chunked batch solve (memory-bounded) ----
    chunk = max(256, floor(8e6 / S));
    for c0 = 1:chunk:M
        cc = c0:min(c0+chunk-1, M);
        [u, v, q, nit] = icgn_batch(ref, Fi, gxf, gyf, x0(cc), y0(cc), ...
            xi, eta, half, du0, dv0, P, Himg, Wimg);
        lin = sub2ind([nY, nX], Ilist(cc), Jlist(cc));
        U(lin) = u; V(lin) = v; Q(lin) = q; NIT(lin) = nit;
    end

    % quality filter (dic_engine never finalizes points below zncc_min)
    bad = Q < P.zncc_min;
    U(bad) = NaN; V(bad) = NaN; NIT(bad) = NaN; Q(bad) = NaN;

    R = struct('U', U, 'V', V, 'Q', Q, 'NIT', NIT, 'gx', gx, 'gy', gy);
end

function [u, v, q, nit] = icgn_batch(ref, Fi, gxf, gyf, x0, y0, ...
    xi, eta, half, du0, dv0, P, Himg, Wimg)
% Vectorized IC-GN for one chunk of n independent points.
    n = numel(x0);
    S = numel(xi);
    u = nan(1, n); v = nan(1, n); q = nan(1, n); nit = nan(1, n);

    % gather reference subsets and their gradients: S x n index matrices
    rows = eta + y0;                             % S x n (implicit expansion)
    cols = xi  + x0;
    lin  = (cols - 1) * Himg + rows;
    f    = ref(lin);
    fz   = f - mean(f, 1);
    fn   = sqrt(sum(fz.^2, 1));                  % 1 x n
    Jx   = gxf(lin);  Jy  = gyf(lin);
    Jxx  = Jx .* xi;  Jxe = Jx .* eta;
    Jyx  = Jy .* xi;  Jye = Jy .* eta;

    % fixed 6x6 Hessian per point (built once, like dic_engine)
    B = {Jx, Jxx, Jxe, Jy, Jyx, Jye};
    H = zeros(6, 6, n);
    for a = 1:6
        for b = a:6
            hab = sum(B{a} .* B{b}, 1);
            H(a, b, :) = hab;
            H(b, a, :) = hab;
        end
    end

    % initial warp: common integer offset (see header note)
    W = repmat([1 0 du0; 0 1 dv0; 0 0 1], 1, 1, n);

    alive = fn > eps;                            % degenerate subsets fail
    act   = find(alive);
    gzL   = zeros(S, n); gnL = zeros(1, n);      % last evaluation per point
    itL   = zeros(1, n);
    for it = 1:P.max_iter
        if isempty(act), break; end
        w11 = reshape(W(1,1,act), 1, []); w12 = reshape(W(1,2,act), 1, []);
        w13 = reshape(W(1,3,act), 1, []); w21 = reshape(W(2,1,act), 1, []);
        w22 = reshape(W(2,2,act), 1, []); w23 = reshape(W(2,3,act), 1, []);
        xw = xi .* w11 + eta .* w12 + (x0(act) + w13);   % S x na
        yw = xi .* w21 + eta .* w22 + (y0(act) + w23);

        % points whose subset left the frame fail (same rule as dic_engine)
        out = min(xw,[],1) < 2 | max(xw,[],1) > Wimg-1 | ...
              min(yw,[],1) < 2 | max(yw,[],1) > Himg-1;
        if any(out)
            alive(act(out)) = false;
            act = act(~out); xw = xw(:, ~out); yw = yw(:, ~out);
            if isempty(act), break; end
        end

        g  = reshape(Fi(yw(:), xw(:)), S, []);
        gz = g - mean(g, 1);
        gn = sqrt(sum(gz.^2, 1));
        dead = gn < eps;
        if any(dead)
            alive(act(dead)) = false;
            keep = ~dead;
            act = act(keep); gz = gz(:, keep); gn = gn(keep);
            if isempty(act), break; end
        end
        gzL(:, act) = gz; gnL(act) = gn; itL(act) = it;

        e = gz .* (fn(act) ./ gn) - fz(:, act);          % ZNSSD residual
        rhs = zeros(6, 1, numel(act));
        for a = 1:6
            rhs(a, 1, :) = sum(B{a}(:, act) .* e, 1);
        end
        dp = pagemldivide(H(:, :, act), rhs);            % 6 x 1 x na
        bad = squeeze(any(~isfinite(dp), 1)).';          % singular Hessian
        if any(bad)
            alive(act(bad)) = false;
            keep = ~bad;
            act = act(keep); dp = dp(:, :, keep);
            if isempty(act), break; end
        end

        % inverse-compositional update W <- W / W(dp), all points at once
        na = numel(act);
        Wdp = zeros(3, 3, na);
        Wdp(1,1,:) = 1 + dp(2,1,:); Wdp(1,2,:) = dp(3,1,:); Wdp(1,3,:) = dp(1,1,:);
        Wdp(2,1,:) = dp(5,1,:);     Wdp(2,2,:) = 1 + dp(6,1,:); Wdp(2,3,:) = dp(4,1,:);
        Wdp(3,3,:) = 1;
        W(:, :, act) = pagemrdivide(W(:, :, act), Wdp);

        % convergence test (identical norm to dic_engine)
        d = reshape(dp, 6, na);
        cnorm = sqrt(d(1,:).^2 + (d(2,:)*half).^2 + (d(3,:)*half).^2 + ...
                     d(4,:).^2 + (d(5,:)*half).^2 + (d(6,:)*half).^2);
        act = act(cnorm >= P.conv_tol);          % converged points retire
    end

    ok = alive;                                  % capped points stay valid,
    u(ok)   = reshape(W(1,3,ok), 1, []);         % exactly like dic_engine
    v(ok)   = reshape(W(2,3,ok), 1, []);
    q(ok)   = sum(fz(:,ok) .* gzL(:,ok), 1) ./ max(fn(ok) .* gnL(ok), eps);
    nit(ok) = itL(ok);
end

function [du, dv] = int_search(ref, cur, y0, x0, half, rng_px)
% Integer ZNCC search around (x0, y0) - the seed only.
% (Verbatim copy of the local function in dic_engine.m; script-local
% functions are not callable across files.)
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
