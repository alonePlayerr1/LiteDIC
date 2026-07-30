function B = calibrate_stereo_ba(cal, opts)
% CALIBRATE_STEREO_BA  Joint stereo bundle adjustment of a two-camera rig.
%
%   B = CALIBRATE_STEREO_BA(cal)          cal from CALIBRATE_CAMERAS
%   B = CALIBRATE_STEREO_BA(cal, opts)
%
%   Refines, in ONE least-squares problem:
%     - both cameras' intrinsics and radial distortion (6 + 6),
%     - the rig's relative pose camera0 -> camera1 (6), shared by all views,
%     - one target pose per view, in camera 0's frame (6 each),
%   minimising the reprojection error in BOTH cameras at once.
%
%   Per-camera Zhang on its own with Zhang's method
%   leaves focal length weakly determined whenever the target was not
%   strongly tilted and not moved in depth: f trades off against object
%   distance, and the reprojection error cannot see the difference. On the
%   ViC3D telephoto rig that produced two cameras 17 % apart in focal length
%   - carrying the same lens - and out-of-plane displacement 4 % low, at a
%   reprojection RMS of 0.14 px. Making the relative pose a SHARED unknown
%   across every view couples the two cameras: the same target points must
%   reproject correctly in both from one rigid baseline, and a focal length
%   that is too long can no longer be absorbed by pushing the target away.
%   See CAL_CONDITIONING for the geometry test that predicts when this
%   matters.
%
%   opts
%     .iters    max LM iterations (default 40)
%     .verbose  print per-iteration RMS (default true)
%     .fixdist  hold both distortion models fixed (default false)
%     .skew     solve for the sensor skew fs (default FALSE - see below)
%     .reject   drop an OBSERVATION whose reprojection residual exceeds
%               median + this many robust sigmas, then refit (default 8;
%               Inf keeps everything). See the note by the rejection.
%
%   Returns B:
%     .K0 .dist0 .K1 .dist1   refined intrinsics / [k1 k2]
%     .R .T                   relative pose camera0 -> camera1 (T in the
%                             units of cal.pitch, i.e. mm for a metric target)
%     .rvec .tvec             per common view, target pose in camera 0
%     .views                  view ids used
%     .rms .rms0 .rms1        reprojection RMS, both cameras and per camera
%     .rms_before             the same, before refinement
%     .dropped                observations rejected as mis-associated
%     .converged              true if LM stopped improving rather than
%                             running out of iterations
%
%   Requires cal.cam(i).views (obj/uv per kept view), which CALIBRATE_CAMERAS
%   stores for exactly this purpose. Base MATLAB, no toolboxes; the Jacobian
%   is assembled sparse, so the per-view pose blocks cost nothing.
%
%   ViC3D rig vs the reference calibration of the same 75 views:
%
%     stage                        fx cam0/cam1   ratio  baseline   RMS
%     per-camera Zhang             17271 / 14780  1.169   405.5 mm  0.14 px
%     + joint adjustment           16339 / 16208  1.008   338.8 mm  0.49 px
%     + rejecting mis-associations 16138 / 16126  1.0007  340.6 mm  0.115 px
%     reference (Vic-3D)           16005 / 15968  1.0023  340.0 mm  0.060 px
%
%   End to end, out-of-plane displacement went from 4.1 % low, to 2.2 % high,
%   to 0.7 % high. Dead ends: the residual split was not missing skew (see
%   REPROJ below); moving the view poses into camera 1's frame does not test
%   anything, it is the same model reparameterised. Actual cause: see the
%   rejection note below.
%
%   Still open: 0.115 px against the reference's 0.060 px. That gap is now
%   symmetric between the cameras, so it looks like plain dot-centre noise
%   rather than a missing model term - a quadratic field over the sensor
%   explains only 0.10 px of it.
%
%   See also CALIBRATE_CAMERAS, CALIBRATE_ZHANG, CAL_CONDITIONING,
%   STEREO_EXTRINSICS, STEREO_GEOMETRY.

    NBASE = 20;   % see nbase() below
    if nargin < 2, opts = struct(); end
    iters   = getf(opts, 'iters', 40);
    verbose = getf(opts, 'verbose', true);
    fixdist = getf(opts, 'fixdist', false);
    skew    = getf(opts, 'skew', false);   % measured to buy nothing here
    reject  = getf(opts, 'reject', 8);     % robust sigmas; Inf = keep all

    if numel(cal.cam) < 2
        error('calibrate_stereo_ba:cams', 'Needs two calibrated cameras.');
    end
    for i = 1:2
        if ~isfield(cal.cam(i), 'views') || isempty(cal.cam(i).views)
            error('calibrate_stereo_ba:views', ...
                ['cal.cam(%d).views is missing - re-run CALIBRATE_CAMERAS ' ...
                 '(it stores the per-view observations the adjustment needs).'], i);
        end
    end

    % ---- views seen by both cameras, in a common order ----
    id0 = [cal.cam(1).views.id];  id1 = [cal.cam(2).views.id];
    [ids, i0, i1] = intersect(id0, id1);
    nv = numel(ids);
    if nv < 3
        error('calibrate_stereo_ba:common', ...
            'Only %d view(s) common to both cameras - need at least 3.', nv);
    end
    % Each camera reprojects the points IT detected, through the SAME view
    % pose - so the two cameras need not have found the same dots, and no
    % point matching between them is required.
    obs = struct('obj', cell(1, nv), 'uv0', [], 'obj1', [], 'uv1', []);
    for k = 1:nv
        obs(k).obj  = cal.cam(1).views(i0(k)).obj;
        obs(k).uv0  = cal.cam(1).views(i0(k)).uv;
        obs(k).obj1 = cal.cam(2).views(i1(k)).obj;
        obs(k).uv1  = cal.cam(2).views(i1(k)).uv;
    end

    % ---- initial parameters ----
    % 7 intrinsics per camera (fx fy cx cy fs k1 k2), 6 for the shared
    % relative pose, 6 per view for the target pose in camera 0's frame
    p = zeros(NBASE + 6*nv, 1);
    p(1:7)   = intr_of(cal.cam(1));
    p(8:14)  = intr_of(cal.cam(2));
    E = stereo_extrinsics(cal);
    p(15:17) = rot2vec(E.R);
    p(18:20) = E.T(:);
    for k = 1:nv
        p(NBASE + (k-1)*6 + (1:3)) = cal.cam(1).rvec{i0(k)}(:);
        p(NBASE + (k-1)*6 + (4:6)) = cal.cam(1).tvec{i0(k)}(:);
    end

    fixed = [];
    if fixdist, fixed = [fixed, 6, 7, 13, 14]; end
    if ~skew,   fixed = [fixed, 5, 12]; end
    free = setdiff(1:numel(p), fixed);

    [ridx, nres] = res_layout(obs);
    r0 = residuals(p, obs, ridx, nres);
    rms_before = sqrt(mean(r0.^2));
    if verbose
        fprintf('stereo BA: %d common views, %d residuals, %d free parameters\n', ...
            nv, nres, numel(free));
        fprintf('  start   RMS %.4f px\n', rms_before);
    end

    % ---- Levenberg-Marquardt ----
    [p, r, converged] = lm(p, obs, ridx, nres, free, iters, verbose);

    % ---- robust rejection of mis-associated points, then refit ----
    % A spurious blob in one calibration view makes INDEX_GRID hand a few
    % dots the wrong target coordinates (gross errors): in view 55 of the
    % ViC3D camera-1 set the median point
    % residual is 0.216 px and the worst is 31.6 px, in the same view, with
    % 141 dots detected where the target has 140. Per-camera Zhang hides it
    % (a free pose per view per camera absorbs a lot); the joint fit cannot,
    % and eight such views dragged camera 1 from 0.16 px to 0.68 px.
    %
    % Rejection is per POINT, not per view: throwing away a whole view would
    % discard ~140 good observations to remove two bad ones. Threshold is loose,
    % only for association failures (orders of magnitude out), not noise.
    B.dropped = 0;
    if isfinite(reject)
        [keepm, ndrop, thr] = point_mask(obs, ridx, r, reject);
        if ndrop > 0 && ndrop < 0.1 * numel(r) / 2
            for k = 1:numel(obs)
                m0 = keepm{k}.c0;  m1 = keepm{k}.c1;
                obs(k).obj  = obs(k).obj(m0, :);   obs(k).uv0 = obs(k).uv0(m0, :);
                obs(k).obj1 = obs(k).obj1(m1, :);  obs(k).uv1 = obs(k).uv1(m1, :);
            end
            B.dropped = ndrop;
            if verbose
                fprintf('  rejecting %d of %d observations over %.3f px\n', ...
                    ndrop, numel(r)/2, thr);
            end
            [ridx, nres] = res_layout(obs);
            % the residual is recomputed from p in the unpack below, over the
            % kept observations, so the one this returns is not needed
            [p, ~, converged] = lm(p, obs, ridx, nres, free, iters, verbose);
        end
    end

    % ---- unpack ----
    B.K0 = k_of(p(1:5));   B.dist0 = p(6:7).';
    B.K1 = k_of(p(8:12));  B.dist1 = p(13:14).';
    B.R  = vec2rot(p(15:17));  B.T = p(18:20);
    B.rvec = cell(1, nv); B.tvec = cell(1, nv);
    for k = 1:nv
        B.rvec{k} = p(NBASE + (k-1)*6 + (1:3));
        B.tvec{k} = p(NBASE + (k-1)*6 + (4:6));
    end
    B.views = ids;
    r = residuals(p, obs, ridx, nres);
    B.rms  = sqrt(mean(r.^2));
    m0 = false(nres, 1);
    for k = 1:nv, m0(ridx{k}.c0) = true; end
    B.rms0 = sqrt(mean(r(m0).^2));
    B.rms1 = sqrt(mean(r(~m0).^2));
    B.rms_before = rms_before;
    B.converged  = converged;
    if verbose
        fprintf(['  done    RMS %.4f px (cam0 %.4f, cam1 %.4f) | was %.4f\n' ...
                 '  fx %.1f / %.1f  (ratio %.4f, was %.4f) | baseline %.3f\n'], ...
            B.rms, B.rms0, B.rms1, rms_before, B.K0(1,1), B.K1(1,1), ...
            B.K0(1,1)/B.K1(1,1), cal.K0(1,1)/cal.K1(1,1), norm(B.T));
    end
end

% ===================== Levenberg-Marquardt =====================

function [p, r, converged] = lm(p, obs, ridx, nres, free, iters, verbose)
    lambda = 1e-3;
    r = residuals(p, obs, ridx, nres);
    prev = sqrt(mean(r.^2));
    converged = false;  cur = prev;
    for it = 1:iters
        J = jacobian(p, obs, ridx, nres, free);
        Jt = J.' * J;  g = J.' * r;
        stepped = false;
        for tries = 1:12
            A = Jt + lambda * spdiags(full(diag(Jt)), 0, numel(free), numel(free));
            dp = -A \ g;
            pn = p; pn(free) = pn(free) + dp;
            rn = residuals(pn, obs, ridx, nres);
            cur = sqrt(mean(rn.^2));
            if cur < prev
                p = pn; r = rn;
                lambda = max(lambda * 0.5, 1e-9);
                stepped = true;
                break;
            end
            lambda = min(lambda * 4, 1e9);
        end
        if verbose && stepped
            fprintf('  it %2d   RMS %.4f px\n', it, cur);
        end
        if ~stepped || abs(prev - cur) < 1e-7
            % stopped because it could not improve, not because it ran out
            % of iterations
            converged = stepped || it > 1;
            break;
        end
        prev = cur;
    end
end

function [keepm, ndrop, thr] = point_mask(obs, ridx, r, k)
% Per-observation residual magnitude, and which ones to keep.
    nv = numel(obs);
    mag = [];
    for q = 1:nv
        mag = [mag; pmag(r, ridx{q}.c0); pmag(r, ridx{q}.c1)]; %#ok<AGROW>
    end
    thr = median(mag) + k * mad_sigma(mag);
    keepm = cell(1, nv);  ndrop = 0;
    for q = 1:nv
        m0 = pmag(r, ridx{q}.c0) <= thr;
        m1 = pmag(r, ridx{q}.c1) <= thr;
        keepm{q} = struct('c0', m0, 'c1', m1);
        ndrop = ndrop + nnz(~m0) + nnz(~m1);
    end
end

function m = pmag(r, idx)
% residuals come stacked [du1 dv1 du2 dv2 ...] per block
    if isempty(idx), m = zeros(0, 1); return; end
    e = reshape(r(idx), 2, []).';
    m = sqrt(sum(e.^2, 2));
end

function s = mad_sigma(x)
    s = 1.4826 * median(abs(x - median(x)));
    if s <= 0, s = std(x); end
    if s <= 0, s = eps; end
end

% ===================== residuals =====================

function [ridx, nres] = res_layout(obs)
    nv = numel(obs); ridx = cell(1, nv); n = 0;
    for k = 1:nv
        a = 2 * size(obs(k).obj, 1);
        b = 2 * size(obs(k).obj1, 1);
        ridx{k}.c0 = n + (1:a);       n = n + a;
        ridx{k}.c1 = n + (1:b);       n = n + b;
        ridx{k}.all = [ridx{k}.c0, ridx{k}.c1];
    end
    nres = n;
end

function r = residuals(p, obs, ridx, nres)
    r = zeros(nres, 1);
    for k = 1:numel(obs)
        [e0, e1] = view_resid(p, obs(k), k);
        r(ridx{k}.c0) = e0;
        r(ridx{k}.c1) = e1;
    end
end

function [e0, e1] = view_resid(p, o, k)
% One view: the target pose is in CAMERA 0's frame, and camera 1 sees the
% same target through the rig's shared relative pose.
    base = nbase() + (k-1)*6;
    rw = p(base+(1:3)); tw = p(base+(4:6));
    Rw = vec2rot(rw);
    e0 = reproj(p(1:7), Rw, tw, o.obj, o.uv0);
    Rr = vec2rot(p(15:17)); Tr = p(18:20);
    e1 = reproj(p(8:14), Rr*Rw, Rr*tw + Tr, o.obj1, o.uv1);
end

function e = reproj(intr, R, t, obj, uv)
% Residuals (predicted - observed), stacked [du1 dv1 du2 dv2 ...]'.
%
% SKEW (intr(5)) is available but OFF by default. This rig's adjustment leaves camera 0 at 0.12 px and camera 1 at
% 0.68 px, and the reference calibration reports fs = -0.98 and +5.44 px -
% which over a normalised half-field of 0.13 is 0.13 and 0.71 px, matching
% those residuals almost exactly. But freeing fs moved the total RMS 0.4899 -> 0.4877 px, left
% camera 1 at 0.68, drove fs to -15.4 px (wrong sign, three times the
% reference) and made the focal-length ratio worse, 1.0081 -> 1.0141. So
% fixed at zero.
    fx = intr(1); fy = intr(2); u0 = intr(3); v0 = intr(4);
    fs = intr(5); k1 = intr(6); k2 = intr(7);
    Pc = R(:,1:2) * obj(:,1:2).' + t;      % planar target, Z = 0
    xn = Pc(1,:) ./ Pc(3,:);
    yn = Pc(2,:) ./ Pc(3,:);
    r2 = xn.^2 + yn.^2;
    d  = 1 + k1*r2 + k2*r2.^2;
    u  = fx * (xn .* d) + fs * (yn .* d) + u0;
    vv = fy * (yn .* d) + v0;
    e  = [u - uv(:,1).'; vv - uv(:,2).'];
    e  = e(:);
end

% ===================== sparse Jacobian =====================

function J = jacobian(p, obs, ridx, nres, free)
% Built from triplets. Camera-0 intrinsics touch only camera-0 residuals,
% camera-1 intrinsics and the relative pose only camera-1 residuals, and a
% view's pose only that view's rows - a dense Jacobian would be 95 % zeros.
    nv = numel(obs);
    col = zeros(1, numel(p));  col(free) = 1:numel(free);

    e0 = cell(1, nv); e1 = cell(1, nv);
    for k = 1:nv, [e0{k}, e1{k}] = view_resid(p, obs(k), k); end

    % Upper bound on the triplet count, filled with a running cursor: growing
    % three arrays inside the parameter loop is what makes a hand-rolled
    % Jacobian slow.
    nnz_max = 0;
    for k = 1:nv, nnz_max = nnz_max + 7*numel(ridx{k}.c0) + 13*numel(ridx{k}.c1); end
    nnz_max = nnz_max + 6 * nres;
    I = zeros(nnz_max, 1); Jc = I; Vv = I;  n = 0;

    for c = 1:nbase()
        if col(c) == 0, continue; end
        h = fdstep(p(c));
        pp = p; pp(c) = pp(c) + h;
        for k = 1:nv
            [f0, f1] = view_resid(pp, obs(k), k);
            if c <= 7
                rows = ridx{k}.c0(:);  vals = (f0 - e0{k}) / h;
            else
                rows = ridx{k}.c1(:);  vals = (f1 - e1{k}) / h;
            end
            q = n + (1:numel(rows));
            I(q) = rows; Jc(q) = col(c); Vv(q) = vals;  n = n + numel(rows);
        end
    end

    for k = 1:nv
        base = nbase() + (k-1)*6;
        rows = [ridx{k}.c0(:); ridx{k}.c1(:)];
        for q6 = 1:6
            c = base + q6;
            if col(c) == 0, continue; end
            h = fdstep(p(c));
            pp = p; pp(c) = pp(c) + h;
            [f0, f1] = view_resid(pp, obs(k), k);
            vals = [(f0 - e0{k})/h; (f1 - e1{k})/h];
            q = n + (1:numel(rows));
            I(q) = rows; Jc(q) = col(c); Vv(q) = vals;  n = n + numel(rows);
        end
    end
    J = sparse(I(1:n), Jc(1:n), Vv(1:n), nres, numel(free));
end

function h = fdstep(x)
    h = 1e-6 * max(1, abs(x));
end

% ===================== small helpers =====================

function n = nbase()
% parameters before the per-view block: 7 + 7 intrinsics + 6 relative pose
    n = 20;
end

function v = intr_of(c)
    v = [c.K(1,1); c.K(2,2); c.K(1,3); c.K(2,3); c.K(1,2); c.dist(1); c.dist(2)];
end

function K = k_of(v)
    K = [v(1) v(5) v(3); 0 v(2) v(4); 0 0 1];
end

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

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
