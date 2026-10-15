function run_stereo_dataset(dataDir, varargin)
% RUN_STEREO_DATASET  Full stereo (3D) DIC on a mixed two-camera folder:
% calibration target views + specimen test pairs in one directory, files
% named <stem>_0.<ext> / <stem>_1.<ext> (camera 0 / camera 1).
%
%   run_stereo_dataset('C:\...\DIC 3D dataset', 'SplitAt', 400)
%
% Pipeline:
%   1. pair discovery: stems that have both _0 and _1 files;
%   2. calibration/test split: stems whose numeric part < 'SplitAt' are
%      calibration views (fast); with 'SplitAt' unset, every stem is
%      classified by actually detecting the dot target (slow but general);
%   3. per-camera Zhang calibration (CALIBRATE_CAMERAS with 'invert' for
%      white-on-black Vic-3D targets) at UNIT pitch + relative pose
%      (STEREO_EXTRINSICS) -> geometry via STEREO_GEOMETRY. NOTE: with an
%      unknown target pitch all 3D output is in PITCH UNITS; recover mm
%      afterwards from a known length or a reference export;
%   4. stereo DIC (STEREO_DIC, unchanged) on the test pairs: first test
%      stem = reference state, 'NStates' further states evenly spaced up
%      to the last;
%   5. summary per state + out-of-plane deflection profile W(s) along the
%      specimen's CENTRE LINE (not an average across its width - a plate
%      bulges in both directions) for comparison with a reference export.
%      W is the component along the surface normal, oriented towards the
%      cameras - the same sign convention as STEREO_VIEW_DATA and the CSV
%      export, so every view of one run reports the same number.
%
% Options (name/value):
%   'SplitAt'   numeric threshold on the stem number separating calibration
%               from test stems (default []: classify by detection)
%   'CalDir'    separate calibration-image folder (default '': cal + test
%               mixed in dataDir). When set, every pair in dataDir is a test
%               pair and the calibration views are read from CalDir - the
%               common layout, e.g. stereo_d_sample\cal_images_7mm.
%   'TestOnly'  true when every pair in dataDir is a specimen state and no
%               target views are present (default []: work it out by
%               detection). Not inferred from a supplied calibration: the
%               detection pass also keeps target views OUT of the test
%               series. Setting it skips a detector pass over the whole
%               series (most of the time on a long run).
%   'Cal'       an ALREADY COMPUTED calibration: a path to a cal_stereo.mat,
%               or the struct CALIBRATE_CAMERAS returns. No target views are
%               needed or looked at and no calibration is solved. Mutually
%               exclusive with 'VicCal' (error if both given).
%   'VicCal'    path to a Vic-3D project.xml: use ITS calibration instead of
%               ours. Validates the DIC + triangulation chain in isolation
%               (same geometry as the reference export) and puts the output
%               directly in mm. Mutually exclusive with 'Cal'.
%   'Invert'    calibration target polarity (default []: auto-detect from the
%               reference cal image). detect_dots wants dark dots on bright
%               paper; a white-on-black target needs invert=true. Set
%               true/false to override the auto-detection.
%   'NStates'   number of deformed states to process besides the reference
%               (default 4, evenly spaced, always including the last;
%               Inf = all)
%   'Pitch'     target pitch if known (default 1 = pitch units)
%   'SubsetSize'/'GridStep'  engine parameters (default 31 / 20)
%   'MaxShift'  clip for the global cross-camera shift, px (default 600)
%   'Roi'       [x y w h] specimen rectangle on the camera-0 reference
%               frame (default: auto-detect the textured strip); the same
%               grid restriction the 2D pipeline uses - without it the
%               grid covers the featureless background around the specimen
%   'Frames'    save a tracking storyboard (default true): one PNG per
%               processed state in litedic_out\frames\ with the tracked
%               cam-0 grid drawn on the image (green = tracked, red x =
%               lost since the reference, yellow = ROI) - the eyeball view
%               of WHERE tracking succeeds/fails (e.g. the weld zone)
%
% Output: printed summary + <dataDir>\litedic_out\stereo_results.mat,
% cal_stereo.mat and PNG previews. Base MATLAB, no toolboxes.
% stereo_results.mat also carries 'qual', the calibration quality gate:
% qual.cond (CAL_CONDITIONING: whether the calibration VIEW GEOMETRY could
% pin the focal length at all - a check no reprojection RMS can make),
% qual.ok (false when the calibration is suspect), qual.checks (cellstr of
% BLOCKING reasons, empty when ok), qual.notes (cellstr of ADVISORY notes
% that do NOT clear ok - e.g. the IAC focal fallback, which was measured
% not to imply a bad calibration), qual.epi_med (median epipolar residual, px) and
% per-camera detection statistics - callers/GUIs should downgrade their
% success status when qual.ok is false (a run can "succeed" on top of a
% garbage calibration and produce meaningless 3D numbers).

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), '6_DIC'));

O.SplitAt    = [];
O.NStates    = 4;
O.Pitch      = 1;
O.SubsetSize = 31;
O.GridStep   = 20;
O.MaxShift   = 600;
O.Roi        = [];
O.Frames     = true; % per-state tracking storyboard PNGs (see header)
O.VicCal     = '';   % path to a Vic-3D project.xml: use ITS calibration
                     % instead of ours (validates the DIC+triangulation
                     % chain in isolation; output directly in mm)
O.Invert     = [];   % calibration target polarity: [] = auto-detect from the
                     % reference cal image (dark dots on bright paper => no
                     % invert; bright dots on dark bg => invert). Set true/false
                     % to force. detect_dots wants black-on-white.
O.CalDir     = '';   % separate calibration-image folder: when set, every
                     % pair in dataDir is a TEST pair and the calibration
                     % views are taken from CalDir (the common layout where
                     % cal images live in their own folder). Default '' =
                     % calibration + test mixed in dataDir (SplitAt/detect).
O.TestOnly   = [];   % true: every pair in dataDir is a test pair and no
                     % target views are looked for (e.g. a folder with one
                     % specimen's frames); skips a detector pass over the
                     % whole series. Default [] = work it out by detection.
O.Cal        = '';   % an ALREADY COMPUTED calibration: a path to a
                     % cal_stereo.mat, or the struct CALIBRATE_CAMERAS
                     % returns. When set, no target views are needed or
                     % looked at and no calibration is solved - this run
                     % uses that geometry as given. See the note below.
for a = 1:2:numel(varargin)
    O.(varargin{a}) = varargin{a+1};
end
% 'Cal' separates correlation errors from calibration errors: same
% calibration on two runs -> any difference is the engine; two calibrations
% on the same frames -> any difference is the calibration. Also the normal
% lab workflow (calibrate once, shoot all day).
if ~isempty(O.Cal) && ~isempty(O.VicCal)
    error('run_stereo_dataset:twocals', ...
        ['''Cal'' and ''VicCal'' both name a calibration to use. Pick one - ' ...
         'silently preferring either would make the output untraceable.']);
end
% A supplied calibration is loaded here first, two decisions below depend
% on what is in it.
calSrc = '';
calGiven = [];
if ~isempty(O.Cal)
    if isstruct(O.Cal)
        calGiven = O.Cal; calSrc = '(in memory)';
    elseif isfile(O.Cal)
        calGiven = load(O.Cal); calSrc = O.Cal;
    else
        error('run_stereo_dataset:calfile', ...
            'No such calibration file: %s', O.Cal);
    end
    for need = {'K0', 'K1', 'dist0', 'dist1'}
        if ~isfield(calGiven, need{1})
            error('run_stereo_dataset:calfields', ...
                ['%s is not a calibration from CALIBRATE_CAMERAS - no ''%s'' ' ...
                 'field.'], calSrc, need{1});
        end
    end
end
haveCal = ~isempty(calGiven) || ~isempty(O.VicCal);

% Output units: a real pitch (mm) or a Vic-3D calibration means mm,
% otherwise pitch units. Computed once here so all printed lines agree.
% A supplied calibration carries its own pitch (the one it was solved at),
% that one wins over this run's 'Pitch' option.
calPitch = O.Pitch;
if ~isempty(calGiven) && isfield(calGiven, 'pitch') && isscalar(calGiven.pitch)
    calPitch = calGiven.pitch;
end
units = tern(calPitch == 1 && isempty(O.VicCal), 'pitch units', 'mm');
outDir = fullfile(dataDir, 'litedic_out');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% -------- 1. discover pairs --------
[stems, nums, exts] = discover_pairs(dataDir);
if isempty(stems)
    error('run_stereo_dataset: no _0/_1 image pairs found in %s', dataDir);
end
fprintf('Found %d synchronized pairs (%s ... %s)\n', numel(stems), ...
    stems{1}, stems{end});

%% -------- 2. calibration / test split --------
if ~isempty(O.CalDir)
    % Calibration images live in a SEPARATE folder: every dataDir pair is a
    % test pair; calibration views come from CalDir.
    cal_caldir = O.CalDir;
    [cal_stems, ~, cal_exts] = discover_pairs(O.CalDir);
    test_stems = stems;  test_exts = exts;
    fprintf('Separate CalDir: %d calibration views in %s\n', numel(cal_stems), O.CalDir);
    fprintf('Test pairs: %d (%s ... %s)\n\n', numel(test_stems), ...
        test_stems{1}, test_stems{end});
else
    cal_caldir = dataDir;
    if ~isempty(O.SplitAt)
        is_cal = nums < O.SplitAt;
    elseif isequal(O.TestOnly, true)
        % Caller says every pair here is a test pair. Not inferred from a
        % supplied calibration: detection also keeps target views out of
        % the test series, a mixed folder still needs that.
        is_cal = false(size(stems));
    else
        fprintf('No SplitAt/CalDir - classifying every stem by target detection (slow)...\n');
        is_cal = false(size(stems));
        for q = 1:numel(stems)
            try
                img = imread(fullfile(dataDir, sprintf('%s_0.%s', stems{q}, exts{q})));
                imd = double(img); imd = max(imd(:)) - imd;  % white-on-black
                D = detect_dots(imd);
                is_cal(q) = D.n >= 30 && size(D.donut_xy, 1) == 3;
            catch
                % unreadable frame: leave it classified as a test pair
            end
        end
    end
    cal_stems  = stems(is_cal);
    cal_exts   = exts(is_cal);
    test_stems = stems(~is_cal);
    test_exts  = exts(~is_cal);
    fprintf('Split: %d calibration views, %d test pairs (%s ... %s)\n\n', ...
        numel(cal_stems), numel(test_stems), test_stems{1}, test_stems{end});
end
if numel(test_stems) < 2
    error('run_stereo_dataset: fewer than 2 test pairs - check SplitAt / CalDir.');
end
if ~haveCal && numel(cal_stems) < 3
    error('run_stereo_dataset: implausible split - check SplitAt / CalDir.');
end

%% -------- 3. calibration + geometry (unit pitch) --------
if ~isempty(O.VicCal)
    G = [];    % built after the ROI/parameters block (needs a matched pair)
    % Vic-3D calibration -> no cal_stereo.mat written. outDir is not wiped
    % between runs, so delete one left by an earlier run, otherwise the app
    % reads it as this run's calibration.
    stale = fullfile(outDir, 'cal_stereo.mat');
    if isfile(stale), delete(stale); end
elseif ~isempty(calGiven)
% Supplied calibration: nothing is solved, no target view read. Geometry
% is rebuilt from the stored intrinsics and poses the same way as after a
% fresh solve.
cal = calGiven;
% Copy it next to the results (the app reads cal_stereo.mat from the
% output folder for the reprojection RMS).
try
    save(fullfile(outDir, 'cal_stereo.mat'), '-struct', 'cal');
catch err
    warning('run_stereo_dataset:calcopy', ...
        'Could not copy the supplied calibration into %s (%s)', outDir, err.message);
end
fprintf('Using the supplied calibration: %s\n', calSrc);
if isfield(cal, 'rms0') && isfield(cal, 'rms1')
    fprintf('  as solved: reprojection RMS %.4f / %.4f px', cal.rms0, cal.rms1);
    if isfield(cal, 'nviews'), fprintf(' on %d views', cal.nviews); end
    fprintf('\n');
end
go = struct();
go.cal     = cal;
go.dicecal = 'THIS_RIG_HAS_NO_DICE_CAL';
G = stereo_geometry(go);
if ~isfinite(G.pitch_mm)
    G.pitch_mm    = 1;
    G.T_mm        = G.T;
    G.baseline_mm = norm(G.T);
end
fprintf(['\nGeometry: baseline %.2f (%s) | toe-in %.2f deg ' ...
    '(spread %.2f) | rig rigidity T_cv %.3f\n\n'], ...
    norm(G.T), units, G.E.ang_deg, G.E.ang_std, G.E.T_cv);
else
co = struct();
co.caldir = cal_caldir;
co.stems  = cal_stems;
% Target polarity: detect_dots wants DARK dots on BRIGHT paper. Invert only
% for white-on-black targets. Auto-detect from the reference cal image
% (invert=true breaks black-on-white targets like stereo_d_sample AMCalB);
% override with the 'Invert' option.
if isempty(O.Invert)
    gcal = imread(fullfile(cal_caldir, sprintf('%s_0.%s', cal_stems{1}, cal_exts{1})));
    gcal = double(gcal); if ndims(gcal) == 3, gcal = mean(gcal, 3); end
    gs = sort(gcal(:));
    lo = gs(max(1, round(0.02 * numel(gs))));
    hi = gs(min(numel(gs), round(0.98 * numel(gs))));
    co.invert = median(gcal(:)) < 0.5 * (lo + hi);   % dark bg => invert
    fprintf('Target polarity auto-detected: invert = %d (%s dots on %s background)\n', ...
        co.invert, tern(co.invert, 'bright', 'dark'), tern(co.invert, 'dark', 'bright'));
else
    co.invert = O.Invert;
end
co.pitch  = O.Pitch;
co.out    = fullfile(outDir, 'cal_stereo.mat');
cal = calibrate_cameras(co);

go = struct();
go.cal     = cal;
go.dicecal = 'THIS_RIG_HAS_NO_DICE_CAL';    % force pitch-unit output
G = stereo_geometry(go);
if ~isfinite(G.pitch_mm)
    % Object points already carry O.Pitch (mm if known, else unit pitch),
    % so T and the triangulation are already in those units - no rescale.
    G.pitch_mm    = 1;
    G.T_mm        = G.T;
    G.baseline_mm = norm(G.T);
end
fprintf(['\nGeometry: baseline %.2f (%s) | toe-in %.2f deg ' ...
    '(spread %.2f) | rig rigidity T_cv %.3f\n\n'], ...
    norm(G.T), units, G.E.ang_deg, G.E.ang_std, G.E.T_cv);
end

%% -------- 4. stereo DIC on the test series --------
ns = min(O.NStates, numel(test_stems) - 1);
% ns deformed states, evenly spaced, plus the reference (as documented for
% 'NStates'). Note: linspace over ns+1 points gives one state too many.
sel = unique(round(linspace(2, numel(test_stems), ns)));
sel = [1, sel(sel > 1)];                     % reference + evenly spaced states
f0 = cell(1, numel(sel)); f1 = cell(1, numel(sel));
for q = 1:numel(sel)
    f0{q} = fullfile(dataDir, sprintf('%s_0.%s', test_stems{sel(q)}, test_exts{sel(q)}));
    f1{q} = fullfile(dataDir, sprintf('%s_1.%s', test_stems{sel(q)}, test_exts{sel(q)}));
end
fprintf('States: %s (reference) + %s\n\n', test_stems{sel(1)}, ...
    strjoin(test_stems(sel(2:end)), ', '));

% Specimen ROI on the camera-0 reference frame: without it the grid and
% seed cover the featureless background around the narrow specimen strip
% (the exact failure mode of the 2D 4pb series). P.roi flows through the
% stereo layer into every dic_engine call unchanged.
I0r = imread(f0{1});
if isinteger(I0r), I0r = double(I0r) / double(intmax(class(I0r)));
else, I0r = double(I0r); end
if ndims(I0r) == 3, I0r = mean(I0r, 3); end
bands = [];
if ~isempty(O.Roi)
    roi = O.Roi;
else
    [roi, bands] = detect_strip(I0r);
end
fprintf('Specimen ROI (cam0 reference): x %d..%d, y %d..%d\n', ...
    roi(1), roi(1)+roi(3), roi(2), roi(2)+roi(4));
% A specimen speckled in separate zones is measured over the bare span
% between them; the nodes there fail and are dropped. Print it, otherwise
% a run with one zone measured looks normal (DIC 3D dataset).
if size(bands, 1) > 1
    txt = arrayfun(@(k) sprintf('%d..%d', bands(k,1), bands(k,2)), ...
        1:size(bands,1), 'UniformOutput', false);
    fprintf(['  %d textured bands merged into it (rows %s) - the bare span ' ...
        'between them\n  is measured over, so expect dead nodes there\n'], ...
        size(bands, 1), strjoin(txt, ', '));
end
fprintf('\n');

P = struct();
P.subset_size     = O.SubsetSize;
P.grid_step       = O.GridStep;
P.search_range    = 40;
P.search_temporal = 40;
P.black_thr       = 0.15;
P.zncc_min        = 0.75;
P.zncc_min_cross  = 0.40;    % 26 deg toe-in lowers attainable cross ZNCC
P.max_iter        = 30;
P.conv_tol        = 1e-3;
P.max_global      = O.MaxShift;
P.roi             = roi;

if ~isempty(O.VicCal)
    I1r = imread(f1{1});
    if isinteger(I1r), I1r = double(I1r) / double(intmax(class(I1r)));
    else, I1r = double(I1r); end
    if ndims(I1r) == 3, I1r = mean(I1r, 3); end
    G = vic_geometry(O.VicCal, I0r, I1r, P);
    fprintf(['Vic-3D calibration in use: baseline %.1f mm | toe-in %.2f deg | ' ...
        'convention %s (epi med %.3f px on the reference pair)\n\n'], ...
        G.baseline_mm, G.E.ang_deg, G.vic.convention, G.vic.epi_ref);
end

t0 = tic;
S = stereo_dic(f0, f1, G, P);
fprintf('stereo_dic: %d states in %.1f s\n\n', numel(sel), toc(t0));

%% -------- 4b. calibration quality gate --------
% A run can complete "successfully" on top of a bad calibration: most target
% views rejected, the IAC focal fallback triggered, epipolar residuals of
% hundreds of px, no error, wrong 3D numbers (the bundled reference rig
% gets 0.1-0.8 px epipolar RMS). Quality checks are collected here, printed
% and saved ('qual' in stereo_results.mat) so the GUI can show a warning
% instead of plain success.
%
% Two levels. BLOCKING ('checks', clears 'ok'): too few indexed views,
% reprojection RMS > 1 px, epipolar residual > 5 px - these came with
% wrong 3D output. ADVISORY ('notes', 'ok' stays true): the solve took an
% unusual path. The IAC focal fallback is advisory: on the DIC 3D dataset
% it fires on both cameras and depth still agrees with Vic-3D to 0.6 %,
% 0.5 um worst-state rms (compare_vic3d); LM refines the focal anyway.
qual = struct('ok', true, 'checks', {{}}, 'notes', {{}}, 'epi_med', NaN, ...
    'epi_worst', NaN, 'epi_worst_state', '', 'nvalid', [], 'cover', [], ...
    'cover_worst', NaN, 'cover_worst_state', '', 'cond', []);
epi_all = arrayfun(@(fr) fr.epi_med, S.frames);
qual.epi_med = median(epi_all(isfinite(epi_all)));
if qual.epi_med > 5   % same bar vic_geometry uses; bundled rig: 0.1-0.8 px
    qual.checks{end+1} = sprintf(['median epipolar residual %.1f px ' ...
        '(a healthy rig stays well under 5 px)'], qual.epi_med);
end
% Median over states hides one bad state (seen: median 0.133 px, one state
% at 56 px, passed). So the worst state is checked separately.
[qual.epi_worst, kworst] = max(epi_all);
if isfinite(qual.epi_worst)
    qual.epi_worst_state = test_stems{sel(kworst)};
    if qual.epi_worst > 5
        qual.checks{end+1} = sprintf(['state %s: epipolar residual %.1f px ' ...
            '(worst of %d states; a healthy pair stays well under 5 px)'], ...
            qual.epi_worst_state, qual.epi_worst, numel(epi_all));
    end
end
% COVERAGE, per state. A state can lose most of its nodes while the
% calibration and the epipolar residuals stay perfect - the images simply
% stopped carrying enough speckle contrast to correlate. sandwich_1_11:
% states at 14057/13793/13663/13822 valid nodes and one at 1892 (light
% dropped mid-series, speckle contrast 23 DN -> 3.7 DN). Launcher only
% checks the last state. With 13% of nodes left, median/peak/profile of
% that state come from a biased subsample.
% Judged against the REFERENCE state, not against the full grid: the
% reference legitimately covers only the speckled part of the ROI, so its
% own count is the right 100%.
nvalid_all = arrayfun(@(fr) fr.nvalid, S.frames);
qual.nvalid = nvalid_all;
qual.cover  = nvalid_all / max(nvalid_all(1), 1);
if numel(nvalid_all) > 1 && nvalid_all(1) > 0
    [qual.cover_worst, kcov] = min(qual.cover(2:end));
    qual.cover_worst_state = test_stems{sel(kcov + 1)};
    % Bands from sandwich_1_11: good states at 92-98% of the reference,
    % so 60% is well outside the spread and 40% is mostly holes.
    if qual.cover_worst < 0.40
        qual.checks{end+1} = sprintf(['state %s: only %.0f%% of the ' ...
            'reference state''s nodes correlated (%d of %d) - that state''s ' ...
            'field is mostly holes and its summary numbers come from a ' ...
            'biased subsample; check the images for lost speckle contrast'], ...
            qual.cover_worst_state, 100*qual.cover_worst, ...
            nvalid_all(kcov + 1), nvalid_all(1));
    elseif qual.cover_worst < 0.70
        qual.notes{end+1} = sprintf(['state %s correlated %.0f%% of the ' ...
            'reference state''s nodes (%d of %d) - noticeably below the ' ...
            'other states'], qual.cover_worst_state, 100*qual.cover_worst, ...
            nvalid_all(kcov + 1), nvalid_all(1));
    end
end
if isfield(S, 'n_shift_fallback') && S.n_shift_fallback > 0
    qual.notes{end+1} = sprintf(['%d state(s) had no credible global ' ...
        'pre-alignment peak and reused the previous accepted shift ' ...
        '(advisory; check those states'' epipolar residuals above)'], ...
        S.n_shift_fallback);
end
if isempty(O.VicCal)
    for ci = 1:numel(cal.cam)
        c = cal.cam(ci);
        rate = c.n_indexed / max(1, c.n_views);
        qual.cam(ci).detect_rate  = rate;
        qual.cam(ci).n_views      = c.n_views;
        qual.cam(ci).n_indexed    = c.n_indexed;
        qual.cam(ci).rms          = c.rms;
        qual.cam(ci).iac_fallback = c.iac_fallback;
        if rate < 0.5
            qual.checks{end+1} = sprintf(['cam%d: only %d/%d calibration ' ...
                'views survived dot detection (%.0f%%)'], ...
                c.cam, c.n_indexed, c.n_views, 100 * rate);
        end
        if c.iac_fallback
            qual.notes{end+1} = sprintf(['cam%d: closed-form intrinsics ' ...
                'failed - centred-principal-point focal fallback used ' ...
                '(advisory on its own; what decides whether the focal ' ...
                'length was recoverable at all is the view geometry, ' ...
                'checked below)'], c.cam);
        end
        if c.rms > 1
            qual.checks{end+1} = sprintf(...
                'cam%d: calibration reprojection RMS %.2f px (> 1 px)', ...
                c.cam, c.rms);
        end
    end
    % View geometry, blocking: a fronto-parallel single-distance target set
    % leaves focal length unobservable and none of the numbers above show
    % it (ViC3D rig: 0.14 px reprojection, out-of-plane still 4% low).
    % See CAL_CONDITIONING.
    qual.cond = cal_conditioning(cal);
    qual.checks = [qual.checks, qual.cond.blocking];
    qual.notes  = [qual.notes,  qual.cond.notes];
    for ci = 1:numel(cal.cam)
        fprintf(['cam%d view geometry: tilt up to %.1f deg | depth spread ' ...
            '%.1f%%%s\n'], cal.cam(ci).cam, qual.cond.tilt_max(ci), ...
            100 * qual.cond.depth_spread(ci), ...
            tern(qual.cond.degenerate(ci), '  <-- DEGENERATE', ''));
    end
end
qual.ok = isempty(qual.checks);
if ~qual.ok
    fprintf('\n*** CALIBRATION QUALITY WARNING - treat 3D results as unreliable ***\n');
    fprintf('  - %s\n', qual.checks{:});
    % Only for these checks; the geometry block prints its own hint and
    % "check the target polarity" would be misleading there.
    if ~isempty(qual.cond) && numel(qual.checks) > numel(qual.cond.blocking)
        fprintf(['Likely causes: too few sharp/complete target views, wrong target\n' ...
                 'polarity (''Invert'' option), wrong dot pitch, or a target style\n' ...
                 'detect_dots cannot handle. Review the calibration image set.\n']);
    end
    fprintf('\n');
end
if ~isempty(qual.notes)
    fprintf('\nCalibration notes (advisory, not a failure):\n');
    fprintf('  - %s\n', qual.notes{:});
    fprintf('\n');
end

%% -------- 5. summary + deflection profiles (orientation-general) --------
% Fit a plane to the reference surface by PCA - exactly as compare_vic3d does
% - so "deflection" is the surface-NORMAL component W = d3D.n and the profile
% runs along the specimen's own principal in-plane axis e1. No assumption
% that the specimen faces the camera or that it bends along image rows.
[nY, nX, ~] = size(S.X_ref_mm);
P3   = reshape(S.X_ref_mm, [], 3);
vfin = all(isfinite(P3), 2);
% Plane is fitted without the outer ring of the valid mask: subsets there
% straddle the mask boundary and triangulate badly (pressurised plate: 615
% nodes, 1.0 %, below -20 mm, worst -72 mm, surface flat to ~0.5 mm RMS).
% One ring of erosion moves the 1st percentile from -20.2 mm to -0.89 mm.
% Erosion width = half a subset in grid nodes, so the whole subset is inside
% the mask (same idea as PLS_STRAIN).
nEdge = max(1, ceil((P.subset_size / 2) / P.grid_step));
vcore = reshape(erode_mask(reshape(vfin, nY, nX), nEdge), [], 1);
fprintf('Reference plane: fitted on %d of %d nodes (%d-node edge band dropped)\n', ...
    nnz(vcore), nnz(vfin), nEdge);
c0   = mean(P3(vcore, :), 1);
[~, ~, Vv] = svd(P3(vcore, :) - c0, 0);
e1  = Vv(:, 1);            % along the beam (largest in-plane extent)
e2  = Vv(:, 2);            % across the beam
nrm = Vv(:, 3);            % surface normal
% SVD normal has arbitrary sign. Normal points towards the cameras, same as
% STEREO_VIEW_DATA and the CSV export, agrees with Vic-3D (otherwise W had
% opposite signs here and in the export, e.g. +5.664 / -5.680 mm, 11-0108).
if nrm(3) > 0, nrm = -nrm; end
fprintf('Reference plane: surface tilt to camera Z = %.1f deg\n', acosd(abs(nrm(3))));
Sc   = reshape((P3 - c0) * e1,  nY, nX);   % along-beam coordinate per node
Tc   = reshape((P3 - c0) * e2,  nY, nX);   % across-beam coordinate per node
Wref = reshape((P3 - c0) * nrm, nY, nX);   % reference out-of-plane residual

prof = struct('stem', {}, 's', {}, 'W', {});
for k = 1:numel(sel)
    D = S.frames(k).D;
    q = D.Q(isfinite(D.Q));
    % nseeds is how many separate fronts the cross field needed. It belongs
    % on this line because it is the only thing that distinguishes "this
    % region is untrackable" from "this region was never reached": both
    % leave NaN, and only the seed count says which.
    % Warn only on one large unreached area, not on the scattered NaNs any
    % fine-grid run has (otherwise it fires on every state).
    ns = '';
    if isfield(D, 'nseeds')
        ns = sprintf(' | seeds %d', D.nseeds);
        if isfield(D, 'unreached_max') && D.unreached_max > 0.05
            ns = sprintf('%s | UNREACHED: one connected %.0f%% of the speckle was never seeded', ...
                ns, 100 * D.unreached_max);
        end
    end
    fprintf('%-10s shift [%4d %4d] | valid %6d | cross ZNCC med %.3f | epi med %.3f px%s\n', ...
        test_stems{sel(k)}, D.shift(1), D.shift(2), S.frames(k).nvalid, ...
        median(q), S.frames(k).epi_med, ns);
    if k >= 2
        d3  = S.frames(k).d3D_mm;
        Wn  = reshape(reshape(d3, [], 3) * nrm, nY, nX);   % normal deflection
        mag = sqrt(sum(d3.^2, 3));
        U0 = S.frames(k).U0; V0 = S.frames(k).V0;
        fprintf(['           temporal med u %+7.2f v %+7.2f px | ' ...
            'med |d3D| %8.4f | mean W(normal) %+8.4f (%s)\n'], ...
            median(U0(isfinite(U0))), median(V0(isfinite(V0))), ...
            median(mag(isfinite(mag))), mean(Wn(isfinite(Wn))), units);
        [sb, wb] = profile_cut(Sc, Tc, Wn);  % centre-line cut along the beam
        p.stem = test_stems{sel(k)};
        p.s = sb;  p.W = wb;
        prof(end+1) = p; %#ok<AGROW>
    end
end
% RMS and 1..99 % band instead of min/max: extremes are a few
% triangulation outliers on the mask edge (pressurised plate: -72 .. +75 mm
% on a surface flat to ~1 mm). Same nodes rlim() keeps off the colorbar.
% Extremes still printed, labelled.
res = Wref(vcore); res = res(isfinite(res));
if isempty(res), res = 0; end
fprintf(['\nReference surface (core, %d-node edge band dropped): ' ...
    'normal residual RMS %.3f (%s), 1..99%% %.3f .. %.3f, ' ...
    'extremes %.3f .. %.3f, %.0f%% of grid valid\n'], nEdge, ...
    sqrt(mean(res.^2)), units, prctile_(res, 1), prctile_(res, 99), ...
    min(res), max(res), 100 * nnz(vfin) / numel(vfin));

%% -------- 5b. tracking storyboard (one PNG per state) --------
% Frame-by-frame eyeball view of the temporal tracking on camera 0: green
% dots = grid nodes tracked from the reference at their CURRENT position,
% red x = nodes valid in the reference but lost by this state (drawn at
% their reference position, since they have no current one). Shows at a
% glance WHERE tracking holds and where it dies (weld glow, blur, no
% speckle). Visual check for the .mat numbers.
if O.Frames
    fdir = fullfile(outDir, 'frames');
    if ~exist(fdir, 'dir'), mkdir(fdir); end
    [Xg, Yg] = meshgrid(S.gx, S.gy);
    vref = all(isfinite(S.X_ref_mm), 3);
    track_overlay(I0r, Xg, Yg, vref, struct('roi', roi, 'title', ...
        sprintf('%s (reference) | %d/%d nodes triangulated', ...
        test_stems{sel(1)}, nnz(vref), numel(vref))), ...
        fullfile(fdir, sprintf('track_1_%s.png', test_stems{sel(1)})));
    for k = 2:numel(sel)
        Ik = imread(f0{k});
        if ndims(Ik) == 3, Ik = mean(double(Ik), 3); end
        U0 = S.frames(k).U0;  V0 = S.frames(k).V0;
        vk = isfinite(U0) & isfinite(V0);
        lost = vref & ~vk;
        track_overlay(Ik, Xg + U0, Yg + V0, vk, struct( ...
            'lost_xy', [Xg(lost), Yg(lost)], 'roi', roi, 'title', ...
            sprintf('%s | tracked %d/%d | lost since ref %d | epi med %.2f px', ...
            test_stems{sel(k)}, nnz(vk), numel(vk), nnz(lost), ...
            S.frames(k).epi_med)), ...
            fullfile(fdir, sprintf('track_%d_%s.png', k, test_stems{sel(k)})));
    end
    fprintf('Tracking storyboard: %s\\track_*.png (one per state)\n', fdir);
end

save(fullfile(outDir, 'stereo_results.mat'), 'S', 'G', 'P', 'O', 'qual', ...
    'cal_stems', 'test_stems', 'sel', 'prof', 'c0', 'e1', 'e2', 'nrm');

f = figure('Visible', 'off', 'Position', [50 50 1300 450]);
subplot(1, 3, 1);
h = imagesc(S.gx, S.gy, Wref); set(h, 'AlphaData', isfinite(Wref));
axis image; colorbar; clim(rlim(Wref));   % robust: a few triangulation
title(sprintf('Reference surface (normal residual, %s)', units));  % outliers
subplot(1, 3, 2);                            % would otherwise crush the scale
Wn_last = reshape(reshape(S.frames(end).d3D_mm, [], 3) * nrm, nY, nX);
h = imagesc(S.gx, S.gy, Wn_last); set(h, 'AlphaData', isfinite(Wn_last));
axis image; colorbar; clim(rlim(Wn_last));
title(sprintf('W = d3D\\cdotn, last state (%s)', units));
subplot(1, 3, 3); hold on;
hp = gobjects(numel(prof), 1);
for q = 1:numel(prof)
    hp(q) = plot(prof(q).s, prof(q).W, 'LineWidth', 1.2, 'DisplayName', prof(q).stem);
end
grid on;
% Legend: MATLAB caps it at 50 entries and with 122 states it covers the
% axes. Past ten states label only first and last (curves are in order).
if numel(prof) <= 10
    legend('Location', 'best', 'Interpreter', 'none');
else
    legend(hp([1 end]), {prof(1).stem, prof(end).stem}, ...
        'Location', 'best', 'Interpreter', 'none');
end
xlabel(sprintf('along-beam coordinate (%s)', units));
ylabel(sprintf('W = d3D\\cdotn (%s)', units));
title('Out-of-plane deflection profiles'); hold off;
print(f, fullfile(outDir, 'stereo_summary.png'), '-dpng', '-r120'); close(f);

% --- field ON the deformed specimen (cam-0, last state) ---
% The imagesc panel above lives on the reference grid rectangle; here the
% same W = d3D.n values are plotted at the TRACKED positions on the last
% state's camera-0 frame, so the color map follows the actual object.
Ik = imread(f0{end});
if ndims(Ik) == 3, Ik = mean(double(Ik), 3); end
[Xg2, Yg2] = meshgrid(S.gx, S.gy);
U0l = S.frames(end).U0;  V0l = S.frames(end).V0;
f = figure('Visible', 'off', 'Position', [50 50 900 900]);
field_overlay(Ik, Xg2 + U0l, Yg2 + V0l, Wn_last, struct('title', ...
    sprintf('W = d3D\\cdotn, %s (%s) - on the deformed shape', ...
    test_stems{sel(end)}, units)));
print(f, fullfile(outDir, 'fields_on_image.png'), '-dpng', '-r150'); close(f);

fprintf('Saved: %s + stereo_summary.png + fields_on_image.png\n', ...
    fullfile(outDir, 'stereo_results.mat'));
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end

function [stems, nums, exts] = discover_pairs(folder)
% Discover synchronized stereo pairs <stem>_0.<ext> / <stem>_1.<ext> in a
% folder, sorted by the trailing number. Any common image extension; macOS
% AppleDouble ('._*') junk filtered. Used for both the test dir and (when
% CalDir is set) a separate calibration dir.
    dd = dir(folder);
    names = {dd(~[dd.isdir]).name};
    stems = {}; nums = []; exts = {};
    for q = 1:numel(names)
        t = regexpi(names{q}, '^(.*)_0\.(tif|tiff|png|jpe?g|bmp)$', 'tokens', 'once');
        if isempty(t), continue; end
        if strncmp(t{1}, '._', 2), continue; end     % macOS AppleDouble junk
        if ~any(strcmpi(names, [t{1} '_1.' t{2}])), continue; end
        stems{end+1} = t{1}; exts{end+1} = t{2};      %#ok<AGROW>
        n = regexp(t{1}, '(\d+)$', 'tokens', 'once');
        if isempty(n), nums(end+1) = NaN; else, nums(end+1) = str2double(n{1}); end %#ok<AGROW>
    end
    % Series first, then the number. A folder holding several series
    % (1-0000.., 11-0000..) used to be sorted by the number alone, which
    % interleaved them: 1-0000, 11-0000, 1-0001, ... One series sorts exactly
    % as before.
    pre = regexprep(stems, '\d+$', '');
    [upre, ~, ip] = unique(pre);
    [~, order] = sortrows([ip(:), nums(:)]);
    stems = stems(order); nums = nums(order); exts = exts(order);
    if numel(upre) > 1
        warning('run_stereo_dataset:series', ['%s holds %d series (%s); they ' ...
            'are taken one after another. Point the run at one series.'], ...
            folder, numel(upre), strjoin(upre, ', '));
    end
end

function [sc, wc] = profile_cut(Sgrid, Tgrid, Wgrid)
% 1-D deflection profile along ONE cut of the specimen - its centre line.
%
% Not a median across the width: a plate bulges across its width too, the
% median reads low (sandwich series: 7.75 mm vs 9.77 mm on the centre line,
% Vic-3D 10.11 mm).
%
% Bin along the beam, at each station fit W locally in both surface
% coordinates and evaluate on the cut. Same as EXPORT_STEREO_DATA and
% SHOW_DEFLECTION, so profile, CSV and plot match.
    sv = Sgrid(:);  tv = Tgrid(:);  wv = Wgrid(:);
    g = isfinite(sv) & isfinite(tv) & isfinite(wv);
    sv = sv(g);  tv = tv(g);  wv = wv(g);
    if isempty(sv), sc = []; wc = []; return; end
    at = median(tv);                       % the centre line across the part
    nb = min(60, max(5, round(numel(sv) / 10)));
    ed = linspace(min(sv), max(sv), nb + 1);
    sc = 0.5 * (ed(1:end-1) + ed(2:end));
    hw = 1.5 * (ed(2) - ed(1));            % overlapping windows: each fit
    wc = nan(1, nb);                       % needs points on both sides
    minc = max(30, round(0.004 * numel(sv)));
    for b = 1:nb
        m = abs(sv - sc(b)) <= hw;
        if nnz(m) < minc, continue; end
        x = sv(m) - sc(b);  u = tv(m) - at;  y = wv(m);
        if min(nnz(x < 0), nnz(x > 0)) < 0.25 * numel(x), continue; end
        A = [ones(numel(x), 1), x, x.^2, u, u.^2, x .* u];
        c = A \ y;
        r = y - A * c;
        s = 1.4826 * median(abs(r - median(r)));
        if s > 0
            k = abs(r - median(r)) <= 3 * s;
            if nnz(k) >= minc, c = A(k, :) \ y(k); end
        end
        wc(b) = c(1);                      % the fit evaluated on the cut
    end
end

function c = rlim(A)
% Robust color limits: the TIGHTER of two rules in each direction - median
% +/- 3 robust-sigma (MAD) and the 1st..99th percentiles of the data.
%
% MAD alone is right for a strain field (a narrow peak with a heavy tail of
% ROI-edge nodes) and wrong for a smooth deflection bowl: there the values
% spread almost uniformly over the range, the MAD is ~0.3 of it, and +/- 3
% sigma goes nearly a full range past the data on both sides (sandwich
% series last state: data [-10.9, 0.3] mm, colourbar [-17.4, +6.4] mm).
% Percentiles keep the limits on the data, MAD still cuts a heavy tail.
% Toolbox-free; degenerate input falls back sensibly.
    a = A(isfinite(A));
    if numel(a) < 5, c = [0 1]; return; end
    m = median(a);
    d = 1.4826 * median(abs(a - m));       % robust sigma from MAD
    if d <= 0, d = std(a); end             % all-equal-ish -> fall back to std
    if d <= 0, c = [m - 1, m + 1]; return; end
    c = [m - 3*d, m + 3*d];
    s = sort(a);  n = numel(s);            % 1st / 99th percentile, toolbox-free
    p = s(max(1, min(n, round([0.01 0.99] * n))));
    c = [max(c(1), p(1)), min(c(2), p(2))];   % the tighter rule in each direction
    if ~(c(2) > c(1)), c = [m - 3*d, m + 3*d]; end   % near-flat data: MAD alone
end

function m = erode_mask(m, k)
% Shrink a logical mask by k nodes in 4-connectivity, toolbox-free. Border
% nodes of the array count as outside, so a mask touching the array edge is
% trimmed there too.
    for it = 1:k
        p = m;
        m(1, :) = false; m(end, :) = false; m(:, 1) = false; m(:, end) = false;
        m(2:end, :)   = m(2:end, :)   & p(1:end-1, :);
        m(1:end-1, :) = m(1:end-1, :) & p(2:end, :);
        m(:, 2:end)   = m(:, 2:end)   & p(:, 1:end-1);
        m(:, 1:end-1) = m(:, 1:end-1) & p(:, 2:end);
    end
end

function v = prctile_(x, q)
% Nearest-rank percentile, toolbox-free (same helper as detect_dots.m).
    x = sort(x(isfinite(x))); n = numel(x);
    if n == 0, v = 0; return; end
    v = x(min(n, max(1, round(q / 100 * n))));
end

function G = vic_geometry(xmlfile, I0, I1, P)
% Build the stereo geometry from a Vic-3D project.xml calibration block.
% The Euler-angle convention of <orientation> is not documented, so every
% plausible convention is assembled and the one with the smallest median
% epipolar residual on an ACTUAL matched reference pair is selected.
% Output units are the units of the Vic translations (mm).
    txt = fileread(xmlfile);
    tk = regexp(txt, ['<camera[^>]*id="(\d)"[^>]*>\s*<intrinsics>([^<]*)' ...
        '</intrinsics>\s*<distortion[^>]*>([^<]*)</distortion>\s*' ...
        '<orientation>([^<]*)</orientation>'], 'tokens');
    if numel(tk) < 2
        error('vic_geometry: could not parse two cameras from %s', xmlfile);
    end
    K = cell(1, 2); dist = cell(1, 2); ang = cell(1, 2); tr = cell(1, 2);
    for i = 1:2
        c = str2double(tk{i}{1}) + 1;
        v = sscanf(tk{i}{2}, '%f');          % [cx cy fx fy skew]
        d = sscanf(tk{i}{3}, '%f');
        o = sscanf(tk{i}{4}, '%f');          % [rx ry rz tx ty tz]
        K{c}    = [v(3) 0 v(1); 0 v(4) v(2); 0 0 1];
        dist{c} = [d(1) d(2)];
        ang{c}  = o(1:3) * pi / 180;
        tr{c}   = o(4:6);
    end
    spc = regexp(txt, 'spc="([\d.]+)"', 'tokens', 'once');
    if ~isempty(spc)
        fprintf('Vic-3D target pitch from project.xml: %s mm\n', spc{1});
    end

    % candidate relative poses (Euler order x world-to-camera vs pose)
    cand = {};
    builders = {@(a) rotz(a(3)) * roty(a(2)) * rotx(a(1)), ...
                @(a) rotx(a(1)) * roty(a(2)) * rotz(a(3))};
    bnames = {'Rz*Ry*Rx', 'Rx*Ry*Rz'};
    for oi = 1:2
        Rc0 = builders{oi}(ang{1}); Rc1 = builders{oi}(ang{2});
        Ra = Rc1 * Rc0.'; Ta = tr{2} - Ra * tr{1};
        cand{end+1} = {Ra, Ta, [bnames{oi} ', world-to-cam']}; %#ok<AGROW>
        R0 = Rc0.'; t0 = -Rc0.' * tr{1};
        R1 = Rc1.'; t1 = -Rc1.' * tr{2};
        Rb = R1 * R0.'; Tb = t1 - Rb * t0;
        cand{end+1} = {Rb, Tb, [bnames{oi} ', cam pose']}; %#ok<AGROW>
    end

    % one real matched pair to arbitrate the convention
    Pc = P;
    if isfield(P, 'zncc_min_cross'), Pc.zncc_min = P.zncc_min_cross; end
    D = stereo_pair_match(I0, I1, Pc);
    [Xg, Yg] = meshgrid(D.gx, D.gy);
    m = isfinite(D.U) & isfinite(D.V);
    uv0 = [Xg(m), Yg(m)];
    uv1 = [Xg(m) + D.U(m), Yg(m) + D.V(m)];

    best = []; bmed = inf; meds = zeros(1, numel(cand));
    for c = 1:numel(cand)
        F = (K{2}.') \ (skew3(cand{c}{2}) * cand{c}{1}) / K{1};
        l = F * [uv0, ones(size(uv0, 1), 1)].';
        e = abs(sum([uv1, ones(size(uv1, 1), 1)].' .* l, 1) ./ ...
            hypot(l(1, :), l(2, :)));
        meds(c) = median(e);
        if meds(c) < bmed, bmed = meds(c); best = c; end
    end
    fprintf('Convention scan (epi med, px): %s\n', ...
        strjoin(arrayfun(@(c) sprintf('%s=%.2f', cand{c}{3}, meds(c)), ...
        1:numel(cand), 'UniformOutput', false), ' | '));
    if bmed > 5
        % Error, not warning: with the right calibration the best convention
        % gives a fraction of a pixel (0.08 and 0.23 px on the two reference
        % rigs), 5 px means wrong geometry. As a warning the run failed
        % minutes later with "seed point did not converge" instead.
        error('vic_geometry:conv', ...
            ['This calibration does not fit these frames: the best of the ' ...
             'four conventions still leaves an epipolar residual of %.2f px ' ...
             'on the reference pair, where a calibration that belongs to ' ...
             'them gives a fraction of a pixel. Two things do this: a ' ...
             'project.xml from a different rig or session, and a reference ' ...
             'pair that is not a specimen frame at all - a calibration-target ' ...
             'view left at the front of the series puts the target, not the ' ...
             'specimen, in front of the cameras.'], bmed);
    end

    G.K0 = K{1}; G.d0 = dist{1};
    G.K1 = K{2}; G.d1 = dist{2};
    G.R = cand{best}{1}; G.T = cand{best}{2}(:);
    G.T_mm = G.T; G.baseline_mm = norm(G.T);
    G.pitch_mm = 1;                          % translations already in mm
    Tx = skew3(G.T);
    G.F = (G.K1.') \ (Tx * G.R) / G.K0;
    G.E = struct('ang_deg', acosd(max(-1, min(1, (trace(G.R) - 1) / 2))), ...
        'ang_std', 0, 'T_cv', 0);
    G.vic = struct('convention', cand{best}{3}, 'epi_ref', bmed, ...
        'epi_all', meds);
end

function S = skew3(t)
    S = [0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0];
end

function R = rotx(a), R = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]; end
function R = roty(a), R = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]; end
function R = rotz(a), R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]; end
