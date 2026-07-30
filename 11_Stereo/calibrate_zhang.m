function C = calibrate_zhang(views, opts)
% CALIBRATE_ZHANG  Single-camera intrinsic calibration by Zhang's method,
% hand-rolled in base MATLAB (no toolboxes).
%
%   C = CALIBRATE_ZHANG(views) calibrates one camera from a set of planar
%   dot-grid views. `views` is a struct array; each element has
%     .obj  [n x 3] planar object points (Z = 0), in pitch units
%     .uv   [n x 2] matching image points, px
%     .id   (optional) an identifier carried into the report
%
%   Pipeline: a normalized-DLT homography per view -> the closed-form
%   intrinsic solution of Zhang (2000) from the image-of-the-absolute-conic
%   constraints -> analytic per-view extrinsics -> a joint Levenberg-
%   Marquardt refinement of the intrinsics (fx, fy, u0, v0), two radial
%   distortion terms (k1, k2) and every view's pose, minimizing the total
%   reprojection error. Views that fail the homography or produce a
%   non-physical pose are dropped and listed.
%
%   Output struct C:
%     C.K        [3 x 3] intrinsic matrix (skew fixed at 0)
%     C.dist     [k1 k2] radial distortion
%     C.rvec/tvec {1 x m} refined per-view rotation vectors / translations
%     C.rms      overall RMS reprojection error, px
%     C.rms_view [m x 1] per-view RMS, px
%     C.kept     indices (into the input views) that were used
%     C.dropped  indices that were dropped
%     C.fx/fy/u0/v0  scalar convenience copies
%     C.iac_fallback true when the closed-form IAC solve failed and the
%                centred-principal-point focal fallback was used instead
%                (a sign of a weakly tilted / degenerate view set)
%
%   See also DETECT_DOTS, INDEX_GRID, CALIBRATE_CAMERAS.

    if nargin < 2, opts = struct(); end
    min_pts = getf(opts, 'min_pts', 20);
    max_it  = getf(opts, 'max_iter', 60);
    h_tol   = getf(opts, 'h_tol',   1);    % max median transfer error, px
    imsize  = getf(opts, 'imsize',  []);   % [H W], enables the focal fallback

    % ---- select usable views and build homographies ----
    kept = []; H = {}; obj = {}; uv = {};
    dropped = [];
    for v = 1:numel(views)
        P = views(v).obj(:, 1:2); q = views(v).uv;
        if size(P,1) < min_pts, dropped(end+1) = v; continue; end %#ok<AGROW>
        Hv = homography(P, q);
        if isempty(Hv) || any(~isfinite(Hv(:))), dropped(end+1) = v; continue; end %#ok<AGROW>
        % point-level cleanup: spurious detections (reflections, clutter)
        % snapped into the lattice corrupt a few correspondences per view;
        % they inflate the joint fit into a uniform residual floor that
        % neither view rejection nor staged LM can undo. Reject homography
        % outlier POINTS and refit, then screen the view as a whole.
        for it = 1:2
            pr = (Hv * [P, ones(size(P,1), 1)].').';
            res = hypot(pr(:,1)./pr(:,3) - q(:,1), pr(:,2)./pr(:,3) - q(:,2));
            good = res <= max(1.0, 4 * median(res));
            if all(good) || nnz(good) < min_pts, break; end
            P = P(good, :); q = q(good, :);
            Hv = homography(P, q);
            if isempty(Hv) || any(~isfinite(Hv(:))), break; end
        end
        if isempty(Hv) || any(~isfinite(Hv(:))) || size(P, 1) < min_pts
            dropped(end+1) = v; continue; %#ok<AGROW>
        end
        % transfer-error screen: a mis-indexed view still yields a formal
        % homography, but one that no longer predicts its own points; it
        % would spoil the global closed-form solution below.
        pr = (Hv * [P, ones(size(P,1), 1)].').';
        pr = pr(:, 1:2) ./ pr(:, 3);
        terr = median(hypot(pr(:,1) - q(:,1), pr(:,2) - q(:,2)));
        if ~isfinite(terr) || terr > h_tol, dropped(end+1) = v; continue; end %#ok<AGROW>
        kept(end+1) = v; %#ok<AGROW>
        H{end+1} = Hv; obj{end+1} = P; uv{end+1} = q; %#ok<AGROW>
    end
    m = numel(kept);
    if m < 3
        error('calibrate_zhang:tooFewViews', ...
            'only %d usable views (need >= 3).', m);
    end

    % ---- closed-form intrinsics from the homographies ----
    Vmat = zeros(2*m, 6);
    for k = 1:m
        Vmat(2*k-1, :) = vij(H{k}, 1, 2);
        Vmat(2*k,   :) = vij(H{k}, 1, 1) - vij(H{k}, 2, 2);
    end
    [~, ~, Vsv] = svd(Vmat, 0);
    b = Vsv(:, end);
    K0 = intrinsics_from_b(b);
    if isempty(K0)
        K0 = intrinsics_from_b(-b);       % null-vector sign ambiguity
    end
    iac_fallback = false;                 % reported in C for quality gating
    if isempty(K0) && ~isempty(imsize)
        % Fallback for weakly tilted (near-degenerate) view sets, where
        % the IAC system is too ill-conditioned for a stable closed form:
        % fix the principal point at the image centre and estimate a
        % common focal length from each homography's orthonormality
        % constraints; the LM refinement below then frees fx, fy, u0, v0.
        K0 = intrinsics_fallback(H, imsize);
        if ~isempty(K0)
            iac_fallback = true;
            % Not a failure: closed-form guess ill-conditioned (weakly tilted
            % views), LM below refines all intrinsics anyway. Printed as a
            % note, not warning(); also in C.iac_fallback.
            fprintf(['calibrate_zhang: closed-form IAC ill-conditioned - using ' ...
                'centred-principal-point focal fallback (f0 = %.0f px); LM ' ...
                'refinement follows.\n'], K0(1,1));
        end
    end
    if isempty(K0)
        error('calibrate_zhang:badB', 'closed-form intrinsics failed.');
    end

    % ---- analytic per-view extrinsics ----
    rvec = cell(1, m); tvec = cell(1, m);
    for k = 1:m
        [R, t] = extrinsics_from_H(K0, H{k});
        rvec{k} = rot2vec(R); tvec{k} = t;
    end

    % ---- LM refinement with post-fit outlier-view rejection ----
    % p = [fx fy u0 v0 k1 k2, (rvec1 tvec1) ... (rvecm tvecm)]
    % A handheld calibration set can contain views whose detection or
    % indexing is subtly wrong: they pass the per-view homography screen
    % yet drag the JOINT fit (runaway distortion, inflated RMS). After
    % each refinement, views whose own RMS exceeds max(1 px, 3x median)
    % are dropped and the fit repeats from the refined parameters.
    intr = [K0(1,1); K0(2,2); K0(1,3); K0(2,3); 0; 0];
    active = 1:m;
    rms_hist = [];
    for round = 1:3
        ma = numel(active);
        p0 = intr;
        for a = 1:ma
            k = active(a);
            p0 = [p0; rvec{k}(:); tvec{k}(:)]; %#ok<AGROW>
        end
        % Staged refinement: first with the distortion FROZEN at zero -
        % on weak-perspective (near-frontal) sets a free k1/k2 can drag
        % the whole fit into a distortion-dominated local minimum (huge
        % k2, inflated focal, displaced principal point). The all-free
        % pass then starts from a geometrically sane point.
        [p, rms_hist] = lm_refine(p0, obj(active), uv(active), ma, max_it, 5:6);
        [p, rh2]      = lm_refine(p,  obj(active), uv(active), ma, max_it);
        rms_hist = [rms_hist, rh2]; %#ok<AGROW> two short passes, joined once
        intr = p(1:6);
        rms_view_a = zeros(ma, 1);
        for a = 1:ma
            base = 6 + (a-1)*6;
            k = active(a);
            rvec{k} = p(base+1:base+3); tvec{k} = p(base+4:base+6);
            e = reproj_resid(intr, rvec{k}, tvec{k}, obj{k}, uv{k});
            rms_view_a(a) = sqrt(mean(sum(reshape(e, 2, []).^2, 1)));
        end
        thr = max(1.0, 3 * median(rms_view_a));
        bad = rms_view_a > thr;
        if ~any(bad) || nnz(~bad) < 3, break; end
        dropped = [dropped, kept(active(bad))]; %#ok<AGROW>
        active = active(~bad);
    end

    % ---- unpack results ----
    fx = intr(1); fy = intr(2); u0 = intr(3); v0 = intr(4);
    k1 = intr(5); k2 = intr(6);
    C.K = [fx 0 u0; 0 fy v0; 0 0 1];
    C.dist = [k1 k2];
    C.fx = fx; C.fy = fy; C.u0 = u0; C.v0 = v0;
    ma = numel(active);
    C.rvec = cell(1, ma); C.tvec = cell(1, ma); C.rms_view = zeros(ma, 1);
    se = 0; sn = 0;
    for a = 1:ma
        k = active(a);
        C.rvec{a} = rvec{k}; C.tvec{a} = tvec{k};
        e = reproj_resid(intr, rvec{k}, tvec{k}, obj{k}, uv{k});
        C.rms_view(a) = sqrt(mean(sum(reshape(e, 2, []).^2, 1)));
        se = se + sum(e.^2); sn = sn + numel(e)/2;
    end
    C.rms = sqrt(se / sn);
    C.kept = kept(active); C.dropped = sort(dropped);
    C.rms_hist = rms_hist;
    C.iac_fallback = iac_fallback;   % true = closed-form IAC failed and the
                                     % centred-principal-point focal fallback
                                     % seeded the fit (weakly tilted view set;
                                     % treat the calibration with suspicion)
end

% ===================== homography (normalized DLT) =====================

function H = homography(XY, uv)
    n = size(XY, 1);
    if n < 4, H = []; return; end
    [X, Tx] = normalize2(XY); [U, Tu] = normalize2(uv);
    A = zeros(2*n, 9);
    for k = 1:n
        x = X(k,1); y = X(k,2); u = U(k,1); w = U(k,2);
        A(2*k-1, :) = [-x -y -1  0 0 0  u*x u*y u];
        A(2*k,   :) = [ 0 0 0  -x -y -1  w*x w*y w];
    end
    [~, ~, V] = svd(A, 0);
    Hn = reshape(V(:, end), 3, 3).';
    H = Tu \ Hn * Tx;
    if abs(H(3,3)) > eps, H = H / H(3,3); end
end

function [Y, T] = normalize2(X)
    c = mean(X, 1); d = X - c;
    s = sqrt(2) / mean(hypot(d(:,1), d(:,2)));
    T = [s 0 -s*c(1); 0 s -s*c(2); 0 0 1];
    Y = (T * [X, ones(size(X,1),1)].').';
    Y = Y(:, 1:2);
end

% ===================== Zhang closed form =====================

function v = vij(H, i, j)
    v = [ H(1,i)*H(1,j); ...
          H(1,i)*H(2,j) + H(2,i)*H(1,j); ...
          H(2,i)*H(2,j); ...
          H(3,i)*H(1,j) + H(1,i)*H(3,j); ...
          H(3,i)*H(2,j) + H(2,i)*H(3,j); ...
          H(3,i)*H(3,j) ].';
end

function K = intrinsics_fallback(H, imsize)
% Focal length from homographies with the principal point pinned at the
% image centre (classic fallback for weak-perspective calibration sets).
% After centring, K = diag(f, f, 1) and B = K^-T K^-1 gives, per view:
%   h1'Bh2 = 0            ->  f^2 = -(h11 h12 + h21 h22) / (h31 h32)
%   h1'Bh1 = h2'Bh2       ->  f^2 = -(|h1p|^2 - |h2p|^2) / (h31^2 - h32^2)
% where hkp are the first two components. The median over all positive
% estimates is robust to the odd degenerate view.
    u0 = imsize(2) / 2; v0 = imsize(1) / 2;
    T = [1 0 -u0; 0 1 -v0; 0 0 1];
    f2 = [];
    for k = 1:numel(H)
        Hc = T * H{k};
        h1 = Hc(:,1); h2 = Hc(:,2);
        d = h1(3) * h2(3);
        if abs(d) > eps
            est = -(h1(1)*h2(1) + h1(2)*h2(2)) / d;
            if isfinite(est) && est > 0, f2(end+1) = est; end %#ok<AGROW>
        end
        d = h1(3)^2 - h2(3)^2;
        if abs(d) > eps
            est = -((h1(1)^2 + h1(2)^2) - (h2(1)^2 + h2(2)^2)) / d;
            if isfinite(est) && est > 0, f2(end+1) = est; end %#ok<AGROW>
        end
    end
    if isempty(f2), K = []; return; end
    f = sqrt(median(f2));
    K = [f 0 u0; 0 f v0; 0 0 1];
end

function K = intrinsics_from_b(b)
    B11 = b(1); B12 = b(2); B22 = b(3); B13 = b(4); B23 = b(5); B33 = b(6);
    den = B11*B22 - B12^2;
    if abs(den) < eps, K = []; return; end
    v0 = (B12*B13 - B11*B23) / den;
    lam = B33 - (B13^2 + v0*(B12*B13 - B11*B23)) / B11;
    if lam <= 0 || B11 <= 0 || lam*B11/den <= 0, K = []; return; end
    alpha = sqrt(lam / B11);
    beta  = sqrt(lam * B11 / den);
    u0 = -B13 * alpha^2 / lam;            % skew set to 0
    K = [alpha 0 u0; 0 beta v0; 0 0 1];
end

function [R, t] = extrinsics_from_H(K, H)
    % K \ H, not inv(K) * h column by column: one triangular solve for all
    % three columns instead of an explicit inverse plus three products -
    % better conditioned, and what the Code Analyzer asks for (MINV).
    KH  = K \ H;
    lam = 1 / norm(KH(:,1));
    r1 = lam * KH(:,1);
    r2 = lam * KH(:,2);
    r3 = cross(r1, r2);
    t  = lam * KH(:,3);
    if t(3) < 0                           % board must be in front
        r1 = -r1; r2 = -r2; r3 = cross(r1, r2); t = -t;
    end
    [U, ~, V] = svd([r1 r2 r3]);          % nearest rotation
    R = U * diag([1 1 sign(det(U*V'))]) * V';
end

% ===================== rotation vector <-> matrix =====================

function R = vec2rot(r)
    th = norm(r);
    if th < 1e-12, R = eye(3); return; end
    k = r / th;
    Kx = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(th)*Kx + (1 - cos(th))*(Kx*Kx);
end

function r = rot2vec(R)
    th = acos(max(-1, min(1, (trace(R) - 1) / 2)));
    if th < 1e-12, r = [0;0;0]; return; end
    r = th / (2*sin(th)) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end

% ===================== projection / residuals =====================

function e = reproj_resid(intr, rv, tv, obj, uv)
% Residuals (predicted - observed), stacked [du1 dv1 du2 dv2 ...]'.
    fx = intr(1); fy = intr(2); u0 = intr(3); v0 = intr(4);
    k1 = intr(5); k2 = intr(6);
    R = vec2rot(rv);
    Pc = R(:,1:2) * obj(:,1:2).' + tv;    % 3 x n  (planar target, Z = 0)
    xn = Pc(1,:) ./ Pc(3,:);
    yn = Pc(2,:) ./ Pc(3,:);
    r2 = xn.^2 + yn.^2;
    d  = 1 + k1*r2 + k2*r2.^2;
    u  = fx * (xn .* d) + u0;
    vv = fy * (yn .* d) + v0;
    e  = [u - uv(:,1).'; vv - uv(:,2).'];
    e  = e(:);
end

% ===================== Levenberg-Marquardt =====================

function [p, rms_hist] = lm_refine(p, obj, uv, m, max_it, fix_idx)
% Joint LM over intrinsics+distortion (6) and per-view poses (6 each).
% The Jacobian is filled block-wise: intrinsic columns are dense, each
% view's pose columns touch only that view's residuals. fix_idx
% (optional) lists parameter indices HELD FIXED during this refinement
% (used to freeze the distortion in the staged solve).
    if nargin < 6, fix_idx = []; end
    npar = numel(p);
    free = setdiff(1:npar, fix_idx);
    nres = 0; ridx = cell(1, m);
    for k = 1:m, ridx{k} = nres + (1:2*size(obj{k},1)); nres = nres + 2*size(obj{k},1); end

    lambda = 1e-3;
    r = residuals(p, obj, uv, m, ridx, nres);
    prev = sqrt(mean(r.^2));
    rms_hist = prev;
    for it = 1:max_it
        J = jacobian(p, obj, uv, m, ridx, nres);
        J = J(:, free);
        Jt = J.' * J; g = J.' * r;
        ok = false;
        for tries = 1:12
            A = Jt + lambda * diag(diag(Jt));
            dp = -A \ g;
            pn = p; pn(free) = pn(free) + dp;
            rn = residuals(pn, obj, uv, m, ridx, nres);
            cur = sqrt(mean(rn.^2));
            if cur < prev
                p = pn; r = rn; prev = cur;
                lambda = max(lambda * 0.5, 1e-9); ok = true; break;
            else
                lambda = min(lambda * 4, 1e9);
            end
        end
        rms_hist(end+1) = prev; %#ok<AGROW>
        if ~ok || (it > 1 && abs(rms_hist(end-1) - prev) < 1e-7), break; end
    end
end

function r = residuals(p, obj, uv, m, ridx, nres)
    r = zeros(nres, 1);
    intr = p(1:6);
    for k = 1:m
        base = 6 + (k-1)*6;
        rv = p(base+1:base+3); tv = p(base+4:base+6);
        r(ridx{k}) = reproj_resid(intr, rv, tv, obj{k}, uv{k});
    end
end

function J = jacobian(p, obj, uv, m, ridx, nres)
    npar = numel(p);
    J = zeros(nres, npar);
    intr = p(1:6);
    r0 = residuals(p, obj, uv, m, ridx, nres);
    % intrinsic columns (dense)
    for c = 1:6
        h = step(p(c));
        pp = p; pp(c) = pp(c) + h;
        rp = residuals(pp, obj, uv, m, ridx, nres);
        J(:, c) = (rp - r0) / h;
    end
    % per-view pose columns (only that view's residuals)
    for k = 1:m
        base = 6 + (k-1)*6;
        rv = p(base+1:base+3); tv = p(base+4:base+6);
        e0 = reproj_resid(intr, rv, tv, obj{k}, uv{k});
        for c = 1:6
            h = step(p(base+c));
            dv = [rv; tv]; dv(c) = dv(c) + h;
            ep = reproj_resid(intr, dv(1:3), dv(4:6), obj{k}, uv{k});
            J(ridx{k}, base+c) = (ep - e0) / h;
        end
    end
end

function h = step(x)
    h = 1e-6 * max(1, abs(x));
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
