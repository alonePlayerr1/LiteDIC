function G = stereo_geometry(opts)
% STEREO_GEOMETRY  Assemble the full stereo-rig geometry for 3D work.
%
%   G = STEREO_GEOMETRY() gathers everything the stereo DIC layer needs
%   into one struct, computing or loading each piece:
%     - per-camera intrinsics + distortion from cal_stereo.mat (recomputed
%       via CALIBRATE_CAMERAS if the file is missing (it is a generated
%       artifact and not tracked),
%     - the relative pose (R, T) camera0 -> camera1 via STEREO_EXTRINSICS
%       (T in pitch units, since calibration used unit pitch),
%     - the physical dot pitch, anchored from DICe's mm-scaled stereo
%       baseline (DICE_STEREO_BASELINE on the tracked stereo_d_spec
%       cal.xml): pitch_mm = baseline_mm(DICe) / baseline_pitch(ours).
%
%   Output struct G:
%     G.K0, G.d0, G.K1, G.d1   intrinsics / [k1 k2] per camera
%     G.R, G.T                 relative pose cam0 -> cam1 (T in pitch units)
%     G.T_mm                   T scaled to millimetres
%     G.pitch_mm               anchored dot pitch, mm
%     G.baseline_mm            |T| in mm (DICe-anchored)
%     G.F                      fundamental matrix (x1' * F * x0 = 0, px)
%     G.E                      full STEREO_EXTRINSICS output (diagnostics)
%
%   G = STEREO_GEOMETRY(opts) accepts:
%     .cal      path to cal_stereo.mat  (default 9_Validation_DICe/cal_stereo.mat)
%               OR an in-memory calibration struct from CALIBRATE_CAMERAS
%               (used as-is; nothing is loaded or saved - this is how the
%               regression test builds its small deterministic geometry)
%     .dicecal  path to DICe cal.xml    (default 9_Validation_DICe/stereo_d_spec/cal.xml)
%
%   See also CALIBRATE_CAMERAS, STEREO_EXTRINSICS, DICE_STEREO_BASELINE,
%   TRIANGULATE_DLT, STEREO_DIC.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    if nargin < 1, opts = struct(); end
    calfile = getf(opts, 'cal', ...
        fullfile(root, '9_Validation_DICe', 'cal_stereo.mat'));
    dicecal = getf(opts, 'dicecal', ...
        fullfile(root, '9_Validation_DICe', 'stereo_d_spec', 'cal.xml'));
    addpath(here);

    % ---- intrinsics (compute once, then cached in the .mat) ----
    if isstruct(calfile)
        S = calfile;                       % in-memory calibration result
    elseif ~isfile(calfile)
        fprintf('stereo_geometry: %s missing - running calibrate_cameras.\n', calfile);
        S = calibrate_cameras(struct('out', calfile));
    else
        S = load(calfile);
    end
    G.K0 = S.K0; G.d0 = S.dist0;
    G.K1 = S.K1; G.d1 = S.dist1;

    % ---- relative pose ----
    % STEREO_EXTRINSICS averages the pose implied by each common view; the
    % bundle adjustment SOLVES for one pose that fits every view at once,
    % together with the intrinsics it has to be consistent with. Prefer it
    % when it ran. E is kept either way - its per-view spread is the rig
    % rigidity diagnostic, which a single fitted pose cannot report.
    E = stereo_extrinsics(S);
    G.E = E;
    if isfield(S, 'ba') && ~isempty(S.ba)
        G.R = S.ba.R; G.T = S.ba.T(:);
        G.pose_from = 'bundle adjustment';
    else
        G.R = E.R; G.T = E.T(:);
        G.pose_from = 'per-view average';
    end

    % ---- physical pitch anchored from DICe's mm baseline ----
    if isfile(dicecal)
        baseline_mm = dice_stereo_baseline(dicecal);
        % against the baseline of the pose actually in use above, not the
        % per-view average, or the anchor would scale a different geometry
        G.pitch_mm    = baseline_mm / norm(G.T);
        G.baseline_mm = baseline_mm;
        G.T_mm        = G.T * G.pitch_mm;
    elseif strcmp(dicecal, 'THIS_RIG_HAS_NO_DICE_CAL')
        % Caller (e.g. run_stereo_dataset) deliberately has no DICe metric
        % anchor: the physical scale comes from the calibration-target pitch
        % instead, applied by the caller. Expected case: print a note, not a
        % warning.
        fprintf(['stereo_geometry: no DICe metric anchor requested; scale ' ...
                 'comes from the target pitch (results in pitch units here, ' ...
                 'converted by the caller).\n']);
        G.pitch_mm = NaN; G.baseline_mm = NaN; G.T_mm = G.T * NaN;
    else
        warning('stereo_geometry:noDiceCal', ...
            '%s not found - metric scale unavailable, using pitch units.', dicecal);
        G.pitch_mm = NaN; G.baseline_mm = NaN; G.T_mm = G.T * NaN;
    end

    % ---- fundamental matrix for epipolar checks (pixel coordinates) ----
    % X_c1 = R*X_c0 + T  =>  x1n' * [T]x*R * x0n = 0
    % =>  F = K1^-T * [T]x * R * K0^-1  with  x1' * F * x0 = 0.
    Tx = [0 -G.T(3) G.T(2); G.T(3) 0 -G.T(1); -G.T(2) G.T(1) 0];
    G.F = (G.K1.') \ (Tx * G.R) / G.K0;
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
