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
%
%   See also DETECT_DOTS, INDEX_GRID, CALIBRATE_CAMERAS.

    if nargin < 2, opts = struct(); end
    min_pts = getf(opts, 'min_pts', 20);
    max_it  = getf(opts, 'max_iter', 60);

    % ---- select usable views and build homographies ----
    kept = []; H = {}; obj = {}; uv = {};
    dropped = [];
    for v = 1:numel(views)
        P = views(v).obj(:, 1:2); q = views(v).uv;
        if size(P,1) < min_pts, dropped(end+1) = v; continue; end %#ok<AGROW>
        Hv = homography(P, q);
        if isempty(Hv) || any(~isfinite(Hv(:))), dropped(end+1) = v; continue; end %#ok<AGROW>
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
    if isempty(K0)
        error('calibrate_zhang:badB', 'closed-form intrinsics failed.');
    end

    % ---- analytic per-view extrinsics ----
    rvec = cell(1, m); tvec = cell(1, m);
    for k = 1:m
        [R, t] = extrinsics_from_H(K0, H{k});
        rvec{k} = rot2vec(R); tvec{k} = t;
    end

    % ---- pack parameters and refine with Levenberg-Marquardt ----
    % p = [fx fy u0 v0 k1 k2, (rvec1 tvec1) ... (rvecm tvecm)]
    p0 = [K0(1,1); K0(2,2); K0(1,3); K0(2,3); 0; 0];
    for k = 1:m, p0 = [p0; rvec{k}(:); tvec{k}(:)]; end %#ok<AGROW>
    [p, rms_hist] = lm_refine(p0, obj, uv, m, max_it);

    % ---- unpack results ----
    fx = p(1); fy = p(2); u0 = p(3); v0 = p(4); k1 = p(5); k2 = p(6);
    C.K = [fx 0 u0; 0 fy v0; 0 0 1];
    C.dist = [k1 k2];
    C.fx = fx; C.fy = fy; C.u0 = u0; C.v0 = v0;
    C.rvec = cell(1, m); C.tvec = cell(1, m); C.rms_view = zeros(m, 1);
    se = 0; sn = 0;
    for k = 1:m
        base = 6 + (k-1)*6;
        rv = p(base+1:base+3); tv = p(base+4:base+6);
        C.rvec{k} = rv; C.tvec{k} = tv;
        e = reproj_resid([fx fy u0 v0 k1 k2]', rv, tv, obj{k}, uv{k});
        C.rms_view(k) = sqrt(mean(sum(reshape(e, 2, []).^2, 1)));
        se = se + sum(e.^2); sn = sn + numel(e)/2;
    end
    C.rms = sqrt(se / sn);
    C.kept = kept; C.dropped = dropped;
    C.rms_hist = rms_hist;
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
    Ki = inv(K);
    h1 = H(:,1); h2 = H(:,2); h3 = H(:,3);
    lam = 1 / norm(Ki * h1);
    r1 = lam * Ki * h1;
    r2 = lam * Ki * h2;
    r3 = cross(r1, r2);
    t  = lam * Ki * h3;
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

function [p, rms_hist] = lm_refine(p, obj, uv, m, max_it)
% Joint LM over intrinsics+distortion (6) and per-view poses (6 each).
% The Jacobian is filled block-wise: intrinsic columns are dense, each
% view's pose columns touch only that view's residuals.
    npar = numel(p);
    nres = 0; ridx = cell(1, m);
    for k = 1:m, ridx{k} = nres + (1:2*size(obj{k},1)); nres = nres + 2*size(obj{k},1); end

    lambda = 1e-3;
    r = residuals(p, obj, uv, m, ridx, nres);
    prev = sqrt(mean(r.^2));
    rms_hist = prev;
    for it = 1:max_it
        J = jacobian(p, obj, uv, m, ridx, nres);
        Jt = J.' * J; g = J.' * r;
        ok = false;
        for tries = 1:12
            A = Jt + lambda * diag(diag(Jt));
            dp = -A \ g;
            pn = p + dp;
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
