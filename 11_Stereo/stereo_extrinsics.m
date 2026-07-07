function E = stereo_extrinsics(cal)
% STEREO_EXTRINSICS  Relative pose (R, T) of camera 1 w.r.t. camera 0 from
% the per-view board poses recovered during single-camera calibration.
%
%   E = STEREO_EXTRINSICS(cal) takes either a cal_stereo.mat path or the
%   result struct from CALIBRATE_CAMERAS directly, and, for every view seen
%   by BOTH cameras, forms the
%   camera-to-camera transform. A synchronized pair sees the same physical
%   board pose, so with world = board:
%       X_c0 = R0*Xb + t0,  X_c1 = R1*Xb + t1
%       => X_c1 = (R1 R0') X_c0 + (t1 - R1 R0' t0)
%   giving R_v = R1 R0' and T_v = t1 - R_v t0 for each view v. A rigid rig
%   makes these consistent across views; they are averaged robustly (the
%   per-component median rotation vector, re-orthonormalized, and the
%   median translation) and the spread is returned as a rigidity check.
%
%   Output struct E:
%     E.R         [3x3] rotation camera0 -> camera1
%     E.T         [3x1] translation camera0 -> camera1 (PITCH units, object
%                 points use unit pitch)
%     E.baseline  |T| (pitch units)
%     E.nviews    number of common views used
%     E.common    the common view numbers
%     E.T_each    [nviews x 3] per-view T (for inspection)
%     E.T_cv      coefficient of variation of |T| across views (rigidity)
%     E.ang_deg   mean relative rotation angle between the cameras (deg)
%     E.ang_std   std of that angle across views (deg)
%
%   See also CALIBRATE_CAMERAS, READ_DICE_CAL, TRIANGULATE_DLT.

    if ischar(cal) || isstring(cal), S = load(cal); else, S = cal; end
    k0 = S.cam(1).kept; k1 = S.cam(2).kept;
    common = intersect(k0, k1);
    n = numel(common);
    if n < 2, error('stereo_extrinsics:tooFew', 'only %d common views.', n); end

    rvecs = zeros(n, 3); Ts = zeros(n, 3); angs = zeros(n, 1);
    for a = 1:n
        v = common(a);
        i0 = find(k0 == v, 1); i1 = find(k1 == v, 1);
        R0 = vec2rot(S.cam(1).rvec{i0}); t0 = S.cam(1).tvec{i0}(:);
        R1 = vec2rot(S.cam(2).rvec{i1}); t1 = S.cam(2).tvec{i1}(:);
        Rv = R1 * R0.';
        Tv = t1 - Rv * t0;
        rvecs(a,:) = rot2vec(Rv).';
        Ts(a,:) = Tv.';
        angs(a) = norm(rot2vec(Rv)) * 180/pi;
    end

    % robust average (median of rotation vectors and translations)
    rbar = median(rvecs, 1).';
    R = vec2rot(rbar);
    [U, ~, V] = svd(R); R = U * diag([1 1 sign(det(U*V'))]) * V';  % clean rotation
    T = median(Ts, 1).';

    E.R = R; E.T = T; E.baseline = norm(T);
    E.nviews = n; E.common = common(:).';
    E.T_each = Ts;
    E.T_cv = std(vecnorm(Ts, 2, 2)) / mean(vecnorm(Ts, 2, 2));
    E.ang_deg = mean(angs); E.ang_std = std(angs);
end

% ---- local rotation helpers (Rodrigues) ----
function R = vec2rot(r)
    th = norm(r);
    if th < 1e-12, R = eye(3); return; end
    k = r(:) / th;
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(th)*K + (1 - cos(th))*(K*K);
end
function r = rot2vec(R)
    th = acos(max(-1, min(1, (trace(R) - 1) / 2)));
    if th < 1e-12, r = [0;0;0]; return; end
    r = th / (2*sin(th)) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end
