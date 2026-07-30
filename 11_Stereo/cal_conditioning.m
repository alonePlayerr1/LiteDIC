function C = cal_conditioning(cal)
% CAL_CONDITIONING  Is this calibration's VIEW GEOMETRY capable of pinning
% the focal length, or was the answer never observable from these images?
%
%   C = CAL_CONDITIONING(cal)   cal from CALIBRATE_CAMERAS
%
%   Zhang's method recovers focal length from the way a plane's projection
%   changes as the plane is TILTED and MOVED IN DEPTH. Shoot every view
%   fronto-parallel at one distance and focal length trades off against
%   object distance exactly: a 10 % longer f with a 10 % further target
%   reprojects onto the same pixels. RMS stays low and f is wrong, so the
%   check has to look at the poses.
%
%   Returns
%     C.tilt_deg{i}     per-view angle between the target normal and
%                       camera i's optical axis, degrees
%     C.tilt_max(i)     largest tilt camera i ever saw
%     C.depth_spread(i) (max Z - min Z) / median Z over the views
%     C.focal_ratio     fx of camera 0 over fx of camera 1
%     C.degenerate(i)   true when camera i's views cannot pin its focal
%     C.blocking        cellstr, empty when the geometry is adequate
%     C.notes           cellstr, advisory
%
%   Thresholds from two measured rigs:
%
%     rig                          tilt max    depth spread   fx0/fx1
%     DICe stereo_d_spec (works)     35.2 deg      9.4 %       1.009
%     ViC3D sandwich (4 % error)     12.5 deg      1.0 %       1.169
%
%   The ViC3D rig is a 16 000 px focal length telephoto pair: the target
%   was never tilted past 12 deg and never moved more than 1 % in depth,
%   and the resulting calibration reads out-of-plane displacement 4 % low
%   while reprojecting at 0.14 px. Its two cameras, which carry the same
%   lens, came out 17 % apart in focal length.
%
%   See also CALIBRATE_CAMERAS, CALIBRATE_ZHANG, RUN_STEREO_DATASET.

    TILT_MIN_DEG = 20;    % below this the plane was never really tilted
    DEPTH_MIN    = 0.05;  % below this the target never moved in depth
    FOCAL_TOL    = 0.10;  % cross-camera focal disagreement worth saying

    n = numel(cal.cam);
    C.tilt_deg = cell(1, n);
    C.tilt_max = nan(1, n);
    C.depth_spread = nan(1, n);
    C.degenerate = false(1, n);
    C.blocking = {};
    C.notes = {};

    for i = 1:n
        c = cal.cam(i);
        nv = numel(c.rvec);
        if nv < 2, continue; end
        t = zeros(1, nv); Z = zeros(1, nv);
        for k = 1:nv
            R = rodrigues(c.rvec{k});
            t(k) = acosd(min(1, abs(R(3,3))));   % target normal vs optical axis
            tv = double(c.tvec{k});
            Z(k) = tv(3);
        end
        C.tilt_deg{i} = t;
        C.tilt_max(i) = max(t);
        mz = median(Z);
        if mz ~= 0, C.depth_spread(i) = (max(Z) - min(Z)) / abs(mz); end

        C.degenerate(i) = C.tilt_max(i) < TILT_MIN_DEG && ...
            C.depth_spread(i) < DEPTH_MIN;
        if C.degenerate(i)
            C.blocking{end+1} = sprintf(['cam%d: the target was never ' ...
                'tilted past %.0f deg and never moved more than %.1f%% in ' ...
                'depth over %d views - focal length is not observable from ' ...
                'this set, so the 3D scale is unreliable even with low ' ...
                'reprojection RMS. Re-shoot with the target tilted ' ...
                '30-45 deg in several directions and at a few distances.'], ...
                c.cam, C.tilt_max(i), 100*C.depth_spread(i), nv);
        end
    end

    C.focal_ratio = NaN;
    if isfield(cal, 'K0') && isfield(cal, 'K1') && cal.K1(1,1) ~= 0
        C.focal_ratio = cal.K0(1,1) / cal.K1(1,1);
        if abs(C.focal_ratio - 1) > FOCAL_TOL
            C.notes{end+1} = sprintf(['the two cameras came out %.0f%% ' ...
                'apart in focal length (fx %.0f vs %.0f) - expected only ' ...
                'if they really carry different lenses'], ...
                100*abs(C.focal_ratio - 1), cal.K0(1,1), cal.K1(1,1));
        end
    end
end

function R = rodrigues(v)
    v = double(v(:)); th = norm(v);
    if th < 1e-12, R = eye(3); return; end
    k = v / th;
    K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * K * K;
end
