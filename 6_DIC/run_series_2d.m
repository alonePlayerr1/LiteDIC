function run_series_2d(dataDir, varargin)
% RUN_SERIES_2D  Incremental 2D DIC over a mono image series with LARGE
% rigid-body motion between frames (hundreds of px per step).
%
%   run_series_2d('C:\data\series1')
%   run_series_2d(dataDir, 'GridStep', 40, 'Roi', [x y w h])
%
% Why this exists: the standard front end (dic_icgn.m) correlates every
% frame against frame 1 with a +/-15 px seed search and a full-frame grid.
% That fails when the specimen occupies a small band of the frame and
% translates by hundreds of px per loading step. This runner generalizes
% the pipeline, per step k-1 -> k:
%   1. rectangular ROI around the specimen, TRACKED from frame 1 by the
%      measured motion (auto-detected once on frame 1, or user-given
%      [x y w h]); a fresh re-detect is used only as a retry fallback;
%   2. coarse integer shift by decimated ZNCC template match, centered on
%      the previous step's shift (temporal prior);
%   3. subpixel fields with the standard engine (dic_engine.m): the coarse
%      shift feeds P.seed_shift (seed search stays +/-15 px) and the seed
%      is placed on the best-textured spot of the ROI (P.seed_xy);
%   4. Lagrangian accumulation: the frame-1 grid points are tracked
%      through every step (step fields interpolated at current point
%      positions), giving TOTAL displacement vs frame 1 and total strains.
%
% Robustness: a step whose field collapses (almost no valid points) is
% retried once with a freshly re-detected ROI; if it still fails, the run
% STOPS, saves everything up to the failure and dumps ROI crops of the
% offending frame pair to litedic_out\ for eyeball diagnosis (motion blur,
% lighting, wrong band, ...). No silent NaN propagation.
%
% Options (name/value):
%   'Ext'         image mask; default '' = auto-detect the common image
%                 extension (jpg/png/tif/tiff/bmp) with the most files;
%                 pass e.g. '*.tif' to force one
%   'Roi'         [x y w h] specimen rectangle, frame-1 px (default: auto)
%   'GridStep'    grid step, px; default 20 (use 40 for a fast first pass)
%   'SubsetSize'  correlation subset, px, odd; default 41
%   'MaxShift'    coarse-search half-range [x y], full-res px; default [200 500]
%   'SelfTest'    run the synthetic large-shift self-test first; default true
%   'Frames'      save a tracking storyboard: one PNG per processed frame
%                 in litedic_out\frames\ with the tracked grid drawn on the
%                 image (green = alive, red x = lost, yellow = ROI). Default
%                 true. This is the eyeball view of the tracking itself -
%                 when a run goes wrong, the storyboard shows WHERE (grid
%                 drifting off, points dying at blur/glare) without digging
%                 into the .mat results.
%   'Live'        show the same overlay in a VISIBLE figure that updates
%                 after every processed frame (default false) - the
%                 real-time view of the tracking. Needs an interactive
%                 MATLAB session (in `matlab -batch` nothing is shown);
%                 combine with 'Frames' freely (each works alone too).
%   'RemoveRBM'   subtract a global affine from the accumulated Lagrangian
%                 field before strain; default false. On a single 2D camera,
%                 out-of-plane tilt/rotation of the specimen appears as an
%                 apparent GLOBAL affine displacement that leaks into strain;
%                 fitting u = a0+a1*x+a2*y (and v likewise) over the surviving
%                 grid and subtracting it isolates the deformation-gradient
%                 (e.g. bending) signal. OFF by default so the raw field and
%                 the paper figures are unchanged; note that it also removes
%                 any genuinely UNIFORM strain (indistinguishable from tilt
%                 with one camera), so it targets gradient-dominated tests.
%
% Output: per-step log printed to the console,
% plus <dataDir>\litedic_out\ with series_2d_results.mat and PNG previews.
% series_2d_results.mat also stores the per-frame accumulated-displacement
% history (cumU_hist, cumV_hist, size [nY nX nStates]; hist_k = the frame
% index of each state, state 1 = reference). That history feeds a virtual
% extensometer as a CURVE over the whole series, e.g.
%   S = load('litedic_out\series_2d_results.mat');
%   E = virtual_extensometer(S.X0(1,:), S.Y0(:,1), S.cumU_hist, ...
%                            S.cumV_hist, [x1 y1 x2 y2]);   % E.opening is 1xnStates
% straddling a crack/weld/bond line reads its opening frame by frame.
% Base MATLAB only - no toolboxes, no GPU.

%% -------- options --------
O.Ext        = '';       % '' = auto-detect (jpg/png/tif/tiff/bmp); or e.g. '*.tif'
O.Roi        = [];
O.GridStep   = 20;
O.SubsetSize = 41;
O.MaxShift   = [200 500];
O.SelfTest   = true;
O.Frames     = true;     % per-frame tracking storyboard PNGs (see header)
O.Live       = false;    % real-time on-screen tracking view (see header)
O.RemoveRBM  = false;    % subtract a global affine before strain (see header)
for a = 1:2:numel(varargin)
    O.(varargin{a}) = varargin{a+1};
end

P.subset_size  = O.SubsetSize;
P.grid_step    = O.GridStep;
P.search_range = 15;      % around the coarse prior - stays small
P.black_thr    = 0.2;
P.zncc_min     = 0.5;
P.max_iter     = 30;
P.conv_tol     = 1e-3;

%% -------- self-test: engine must recover a LARGE known shift --------
if O.SelfTest
    fprintf('=== Self-test: rigid shift (+60.30, -40.20) px with coarse prior (60, -40) ===\n');
    rng(0);
    img = conv2(rand(420), gk(2.0), 'same');
    img = img(11:410, 11:410);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 60.30; TV = -40.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    sh = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    % The self-test verifies the ENGINE (sign/axis/interpolation), so its
    % subset/step are PINNED to values matched to its own synthetic pattern
    % (speckle sigma 2.0 -> subset 41, as in RUN_TESTS). They must NOT
    % inherit the caller's SubsetSize: a subset tuned to fine real-world
    % speckle (e.g. 15 px for a ~4 px speckle) is too small for THIS coarse
    % pattern and fails the test spuriously and blocks a valid run.
    Pt = P; Pt.subset_size = 41; Pt.grid_step = 10; Pt.seed_shift = [60 -40];
    Rt = dic_engine(img, sh, Pt);
    eu = mean(Rt.U(:), 'omitnan') - TU;  ev = mean(Rt.V(:), 'omitnan') - TV;
    fprintf('  error u: %+.5f px, error v: %+.5f px, valid pts: %d\n', ...
        eu, ev, nnz(~isnan(Rt.U)));
    if abs(eu) < 0.05 && abs(ev) < 0.05
        fprintf('  SELF-TEST PASSED\n\n');
    else
        error('run_series_2d: self-test FAILED - do not trust the results below.');
    end
end

%% -------- load series --------
names = list_series(dataDir, O.Ext);
if isempty(names)
    error('No image files (%s) in %s', tern(isempty(O.Ext), 'auto', O.Ext), dataDir);
end
nF = numel(names);
fprintf('Series: %s  (%d frames: %s ... %s)\n', dataDir, nF, names{1}, names{end});
rd = @(f) read_gray(fullfile(dataDir, f));

prev = rd(names{1});
[H, W] = size(prev);
fprintf('Frame size: %d x %d px | grid step %d px, subset %d px\n', ...
    W, H, P.grid_step, P.subset_size);

outDir = fullfile(dataDir, 'litedic_out');
if ~exist(outDir, 'dir'), mkdir(outDir); end
framesDir = fullfile(outDir, 'frames');
if O.Frames && ~exist(framesDir, 'dir'), mkdir(framesDir); end
hLive = gobjects(0);      % live-view figure, created lazily (see overlay calls)

if ~isempty(O.Roi)
    roi1 = O.Roi;
    fprintf('Manual ROI (frame 1): [%d %d %d %d]\n\n', roi1);
else
    roi1 = detect_band(prev, W, H);
    fprintf(['Auto-ROI (frame 1): cols %d..%d (%.0f%%W), rows %d..%d (%.0f%%H), ' ...
        'tracked by measured motion\n\n'], ...
        roi1(1), roi1(1)+roi1(3), 100*roi1(3)/W, ...
        roi1(2), roi1(2)+roi1(4), 100*roi1(4)/H);
end

ds = max(1, ceil(max(H, W) / 1000));       % decimation for the coarse search
sRc = ceil(O.MaxShift(1) / ds);            % coarse search, cols (X)
sRr = ceil(O.MaxShift(2) / ds);            % coarse search, rows (Y)

%% -------- incremental loop --------
first = rd(names{1});     % kept only for the final quiver plot
trk_x = 0; trk_y = 0;     % measured cumulative motion of the ROI, px
prior = [0 0];            % previous step's coarse shift (temporal prior)
stats = struct('coarse', {}, 'coarse_q', {}, 'roi_y', {}, 'npts', {}, ...
    'zncc', {}, 'step_u', {}, 'step_v', {}, 'alive', {}, ...
    'cum_u', {}, 'cum_v', {}, 't', {});
X0 = []; % frame-1 grid (defined after the first step)
aborted = false;
fprintf('%-4s %-14s %12s %6s %6s %7s %9s %9s %7s %10s %10s %6s\n', 'k', 'frame', ...
    'coarse[x,y]', 'cQ', 'roiY', 'pts', 'step_u', 'step_v', 'ZNCC', ...
    'cum_u', 'cum_v', 'sec');

for k = 2:nF
    t0 = tic;
    cur = rd(names{k});
    if ~isequal(size(cur), [H, W])   % a stray thumbnail/crop in the folder
        error(['run_series_2d: frame %s is %dx%d px but the series is ' ...
            '%dx%d - remove mixed-size files or pass a stricter ''Ext''.'], ...
            names{k}, size(cur, 2), size(cur, 1), W, H);
    end

    % --- (1) ROI in the prev frame = frame-1 ROI tracked by measured motion
    roi = clamp_roi([roi1(1) + round(trk_x), roi1(2) + round(trk_y), ...
                     roi1(3), roi1(4)], W, H);

    % --- (2)+(3) coarse shift + subpixel step field, with one retry ---
    [R, sx, sy, cq, ok] = solve_step(prev, cur, roi, P, ds, sRr, sRc, prior);
    if ~ok && isempty(O.Roi)
        roi_alt = detect_band(prev, W, H);
        fprintf('  step %d: collapsed on tracked ROI (rows %d..%d) - retry with re-detected rows %d..%d\n', ...
            k, roi(2), roi(2)+roi(4), roi_alt(2), roi_alt(2)+roi_alt(4));
        roi = roi_alt;
        [R, sx, sy, cq, ok] = solve_step(prev, cur, roi, P, ds, sRr, sRc, prior);
    end
    if ~ok
        fprintf(['\nSTEP %d FAILED (%s -> %s): field collapsed twice ' ...
            '(coarse [%d,%d], ZNCC %.2f).\n'], k, names{k-1}, names{k}, sx, sy, cq);
        dump_pair(prev, cur, roi, [sx sy], outDir, k);
        fprintf('ROI crops of the failing pair saved to %s\\debug_step%02d_*.png - inspect for blur/lighting.\n', ...
            outDir, k);
        aborted = true;
        break;
    end

    % --- (4) Lagrangian accumulation on the frame-1 grid ---
    if isempty(X0)
        [X0, Y0] = meshgrid(R.gx, R.gy);
        Xc = X0; Yc = Y0;
        cumU = zeros(size(X0)); cumV = zeros(size(X0));
        alive = true(size(X0));
        % per-frame accumulated-displacement history on the frame-1 grid, so
        % a virtual extensometer (virtual_extensometer.m) can be read as a
        % CURVE over the series, not just on the final frame. Snapshot-only:
        % it copies cumU/cumV and changes no computation. State 1 is the
        % reference (zero displacement); one state is appended per step below.
        cumU_hist = zeros([size(X0), 0]);  cumV_hist = cumU_hist;
        cumU_hist(:,:,end+1) = cumU;  cumV_hist(:,:,end+1) = cumV; %#ok<AGROW>
        hist_k = 1;                          % frame index of each history state
        if O.Frames || O.Live    % storyboard frame 1: the initial grid
            inf1 = struct('roi', roi1, 'title', sprintf( ...
                'frame 1  %s | initial grid: %d pts', names{1}, numel(X0)));
            if O.Live
                if isempty(hLive) || ~ishghandle(hLive)
                    hLive = figure('Visible', 'on', 'NumberTitle', 'off', ...
                        'Name', 'LiteDIC live tracking');
                end
                inf1.fig = hLive;
            end
            track_overlay(first, X0, Y0, true(size(X0)), inf1, ...
                tern(O.Frames, fullfile(framesDir, 'track_001.png'), ''));
        end
    end
    u = interp2(R.gx, R.gy, R.U, Xc, Yc, 'linear', NaN);
    v = interp2(R.gx, R.gy, R.V, Xc, Yc, 'linear', NaN);
    died = alive & (~isfinite(u) | ~isfinite(v));
    cumU(died) = NaN;  cumV(died) = NaN;
    alive = alive & ~died;
    Xc(alive) = Xc(alive) + u(alive);
    Yc(alive) = Yc(alive) + v(alive);
    cumU(alive) = cumU(alive) + u(alive);
    cumV(alive) = cumV(alive) + v(alive);
    % Appended, not preallocated: number of states is only known after the
    % loop. Cheap compared to the DIC solve.
    cumU_hist(:,:,end+1) = cumU;  cumV_hist(:,:,end+1) = cumV; %#ok<AGROW>
    hist_k(end+1) = k; %#ok<AGROW>

    % --- track the ROI + temporal prior for the next step ---
    mu = mean(R.U(:), 'omitnan');  mv = mean(R.V(:), 'omitnan');
    trk_x = trk_x + mu;  trk_y = trk_y + mv;
    prior = [sx sy];

    s.coarse = [sx sy];  s.coarse_q = cq;  s.roi_y = roi(2);
    s.npts   = nnz(~isnan(R.U));
    s.zncc   = mean(R.Q(:), 'omitnan');
    s.step_u = mean(u(alive));  s.step_v = mean(v(alive));
    s.alive  = nnz(alive);
    s.cum_u  = mean(cumU(alive));  s.cum_v = mean(cumV(alive));
    s.t      = toc(t0);
    stats(end+1) = s; %#ok<AGROW>
    fprintf('%-4d %-14s [%4d,%5d] %6.2f %6d %7d %9.2f %9.2f %7.3f %10.1f %10.1f %6.1f\n', ...
        k, names{k}, sx, sy, cq, roi(2), s.npts, s.step_u, s.step_v, ...
        s.zncc, s.cum_u, s.cum_v, s.t);

    if O.Frames || O.Live    % storyboard: tracked grid at its CURRENT
                             % positions, dead points (red x) frozen where
                             % they were last seen; optional live view
        infk = struct('lost', ~alive, ...
            'roi', [roi(1) + mu, roi(2) + mv, roi(3), roi(4)], ...
            'title', sprintf(['frame %d  %s | coarse [%d,%d] | ZNCC %.3f | ' ...
            'alive %d/%d | cum [%.1f, %.1f] px'], k, names{k}, sx, sy, ...
            s.zncc, s.alive, numel(alive), s.cum_u, s.cum_v));
        if O.Live
            if isempty(hLive) || ~ishghandle(hLive)
                hLive = figure('Visible', 'on', 'NumberTitle', 'off', ...
                    'Name', 'LiteDIC live tracking');
            end
            infk.fig = hLive;
        end
        track_overlay(cur, Xc, Yc, alive, infk, ...
            tern(O.Frames, fullfile(framesDir, sprintf('track_%03d.png', k)), ''));
    end

    prev = cur;
end

if isempty(X0)
    error('run_series_2d: no step completed - nothing to save.');
end

%% -------- optional rigid-body-motion (global-affine) removal --------
% Fit u = a0+a1*X+a2*Y and v = b0+b1*X+b2*Y over the surviving frame-1 grid
% and subtract. This removes the apparent global affine that out-of-plane
% tilt/rotation leaks into a single-camera field, leaving the deformation-
% gradient (bending) signal. OFF by default (raw field unchanged); it also
% removes any uniform strain, which one camera cannot separate from tilt.
if O.RemoveRBM
    g = alive & isfinite(cumU) & isfinite(cumV);
    if nnz(g) < 6
        warning('run_series_2d:rbm', ...
            'RemoveRBM: too few valid points (%d) to fit an affine - skipped.', nnz(g));
    else
        A  = [ones(nnz(g), 1), X0(g), Y0(g)];
        au = A \ cumU(g);
        av = A \ cumV(g);
        cumU(g) = cumU(g) - A * au;
        cumV(g) = cumV(g) - A * av;
        fprintf(['RBM removal ON: global affine subtracted ' ...
            '(u: a0 %+.3g, du/dx %+.3g, du/dy %+.3g | ' ...
            'v: b0 %+.3g, dv/dx %+.3g, dv/dy %+.3g)\n'], ...
            au(1), au(2), au(3), av(1), av(2), av(3));
    end
end

%% -------- total strains (PLS on the accumulated Lagrangian field) --------
% Green-Lagrange (use_green = true): this driver targets LARGE-deflection
% series where the specimen rotates a lot (e.g. 4pb, thousands of px). Linear
% strain would read that rigid-body rotation as apparent strain; the
% Green-Lagrange measure is rotation-invariant, so the field shows real
% material strain. (The small-strain plate path in dic_icgn.m / the paper
% figure stays linear - unchanged.)
strain_win = max(50, 3 * P.grid_step);
[exx, eyy, exy] = pls_strain(cumU, cumV, P.grid_step, strain_win, true);
fprintf(['\nTotal strain (Green-Lagrange PLS, win %d px): ' ...
    'exx median %+.2e, eyy median %+.2e\n'], ...
    strain_win, median(exx(:), 'omitnan'), median(eyy(:), 'omitnan'));
fprintf('Survived %d / %d grid points%s.\n', nnz(alive), numel(alive), ...
    tern(aborted, ' (run ABORTED early - partial results)', ''));

%% -------- save --------
save(fullfile(outDir, 'series_2d_results.mat'), 'X0', 'Y0', 'Xc', 'Yc', ...
    'cumU', 'cumV', 'cumU_hist', 'cumV_hist', 'hist_k', 'alive', ...
    'exx', 'eyy', 'exy', 'stats', 'P', 'O', 'names', 'aborted');

% --- orientation-free scalar fields (no vertical-plate assumption) ---
% |displacement| and the strain INVARIANTS transform identically under any
% specimen mounting (see test_flip_invariance.m), so they are what the
% headline panels show. Raw components (exx,eyy,exy, u,v) are kept as a
% secondary figure for anyone who needs a specific axis.
dispmag = hypot(cumU, cumV);                        % total displacement, px
pa = (exx + eyy) / 2;
pb = sqrt(((exx - eyy) / 2).^2 + exy.^2);
e1 = pa + pb;  e2 = pa - pb;                         % principal strains
vm = sqrt(e1.^2 - e1.*e2 + e2.^2);                  % von Mises equiv. strain

% Robust color limits (2..98 percentile): large rigid-body rotation in a
% big-deflection series leaves a few extreme strain outliers at the ROI
% edge that otherwise swamp the colorbar and hide the real field.
f = figure('Visible', 'off', 'Position', [50 50 1400 500]);
subplot(1,3,1); h = imagesc(X0(1,:), Y0(:,1), dispmag);
set(h, 'AlphaData', isfinite(dispmag)); axis image; colorbar; clim(rlim(dispmag));
title('|displacement| total, px');
subplot(1,3,2); h = imagesc(X0(1,:), Y0(:,1), e1);
set(h, 'AlphaData', isfinite(e1)); axis image; colorbar; clim(rlim(e1));
title('Max principal strain \epsilon_1 (PLS)');
subplot(1,3,3); h = imagesc(X0(1,:), Y0(:,1), vm);
set(h, 'AlphaData', isfinite(vm)); axis image; colorbar; clim(rlim(vm));
title('von Mises equiv. strain');
print(f, fullfile(outDir, 'fields_total.png'), '-dpng', '-r120'); close(f);

% secondary: raw tensor components (orientation-dependent, for reference)
f = figure('Visible', 'off', 'Position', [50 50 1400 500]);
subplot(1,3,1); h = imagesc(X0(1,:), Y0(:,1), exx);
set(h, 'AlphaData', isfinite(exx)); axis image; colorbar; clim(rlim(exx));
title('\epsilon_{xx} (PLS)');
subplot(1,3,2); h = imagesc(X0(1,:), Y0(:,1), eyy);
set(h, 'AlphaData', isfinite(eyy)); axis image; colorbar; clim(rlim(eyy));
title('\epsilon_{yy} (PLS)');
subplot(1,3,3); h = imagesc(X0(1,:), Y0(:,1), exy);
set(h, 'AlphaData', isfinite(exy)); axis image; colorbar; clim(rlim(exy));
title('\epsilon_{xy} (PLS)');
print(f, fullfile(outDir, 'fields_components.png'), '-dpng', '-r120'); close(f);

% --- fields ON the deformed specimen (last processed frame) ---
% The panels above live on the frame-1 grid, which is always an
% axis-aligned rectangle - a bent beam still shows as a straight band.
% This figure plots the SAME fields at the tracked current positions
% (Xc, Yc) over the last frame, so the color map follows the actual
% deformed shape of the object.
f = figure('Visible', 'off', 'Position', [50 50 1500 600]);
subplot(1, 2, 1);   % prev = last SUCCESSFULLY processed frame (on an
                    % aborted run cur would be the failing one)
field_overlay(prev, Xc, Yc, dispmag, struct('title', ...
    '|displacement| total, px - on the deformed shape'));
subplot(1, 2, 2);
field_overlay(prev, Xc, Yc, vm, struct('title', ...
    'von Mises equiv. strain - on the deformed shape'));
print(f, fullfile(outDir, 'fields_on_image.png'), '-dpng', '-r150'); close(f);

f = figure('Visible', 'off');
imagesc(first); colormap(gray); axis image; hold on;
st = max(1, round(size(X0, 2) / 40));
ii = 1:2:size(X0,1); jj = 1:st:size(X0,2);
quiver(X0(ii,jj), Y0(ii,jj), cumU(ii,jj), cumV(ii,jj), 1, 'r');
title('Total displacement on frame 1'); hold off;
print(f, fullfile(outDir, 'quiver_total.png'), '-dpng', '-r120'); close(f);

f = figure('Visible', 'off');
plot(2:(1+numel(stats)), hypot([stats.cum_u], [stats.cum_v]), 'o-', 'LineWidth', 1.5);
grid on; xlabel('frame'); ylabel('|mean total displacement|, px');
title('Displacement growth over the series');
print(f, fullfile(outDir, 'evolution.png'), '-dpng', '-r120'); close(f);

fprintf('Saved: %s\n', fullfile(outDir, 'series_2d_results.mat'));
fprintf(['PNG previews: fields_total.png (invariants), fields_components.png, ' ...
    'fields_on_image.png (deformed shape), quiver_total.png, evolution.png\n']);
if O.Frames
    fprintf('Tracking storyboard: %s\\track_*.png (one per frame)\n', framesDir);
end
end

%% ==================== local functions ====================
function [R, sx, sy, cq, ok] = solve_step(prev, cur, roi, P, ds, sRr, sRc, prior)
% One incremental step: coarse shift + engine solve on the given ROI.
% ok = false when the field collapses (seed landed on junk / blur / ...).
    [H, W] = size(prev);
    [sx, sy, cq] = coarse_shift(prev, cur, roi, ds, sRr, sRc, prior);
    Pk = P;
    Pk.roi        = roi;
    Pk.seed_shift = [sx sy];
    Pk.seed_xy    = best_texture_xy(prev, roi, W, H);
    R = [];
    try
        R = dic_engine(prev, cur, Pk);
    catch
        ok = false;  return;
    end
    ok = nnz(~isnan(R.U)) >= max(20, 0.01 * numel(R.U));
end

function c = rlim(A)
% Robust color limits: median +/- 3 robust-sigma (MAD-based). This is
% resistant to the HEAVY tail of ROI-edge / low-coverage strain nodes that
% a plain percentile still lets through (Green-Lagrange squares gradients,
% so edge noise is large), so the colorbar zooms onto the real field.
% Toolbox-free. Falls back sensibly on degenerate input.
    a = A(isfinite(A));
    if numel(a) < 5, c = [0 1]; return; end
    m = median(a);
    d = 1.4826 * median(abs(a - m));       % robust sigma from MAD
    if d <= 0, d = std(a); end             % all-equal-ish -> fall back to std
    if d <= 0, c = [m - 1, m + 1]; return; end
    c = [m - 3*d, m + 3*d];
end

function names = list_series(dataDir, ext)
% List image frames in dataDir, sorted, with dotfiles / macOS AppleDouble
% ('._*') junk removed. If ext is given ('*.JPG', '.tif', 'tif', ...) it is
% used verbatim; otherwise the common image extension with the most files
% is auto-selected, so a JPG, TIFF or PNG series all "just work".
    if ~isempty(ext)
        exts = { regexprep(lower(ext), '^\*?\.?', '') };  % '*.JPG' -> 'jpg'
    else
        exts = {'jpg','jpeg','png','tif','tiff','bmp'};
    end
    best = {}; bestn = -1;
    for k = 1:numel(exts)
        d  = dir(fullfile(dataDir, ['*.' exts{k}]));
        nm = {d.name};
        nm = nm(~startsWith(nm, '.'));                    % drop dotfiles/junk
        if numel(nm) > bestn, bestn = numel(nm); best = nm; end
    end
    names = sort_frames(best);
end

function names = sort_frames(names)
% Chronological sort of frame names. Real cameras/DAQs often number frames
% WITHOUT zero-padding (img_2.jpg ... img_10.jpg): a plain lexicographic
% sort puts img_10 before img_2 and the series is processed out of order.
% When every name carries a trailing number before the extension, sort by
% that number (MATLAB sort is stable, so equal numbers keep dir order);
% otherwise fall back to the lexicographic sort. Zero-padded series sort
% identically either way, so existing datasets are unaffected.
    if isempty(names), return; end
    t = regexp(names, '(\d+)\.[^.]*$', 'tokens', 'once');
    has = ~cellfun(@isempty, t);
    if all(has)
        num = cellfun(@(c) str2double(c{1}), t);
        [~, o] = sort(num);
        names = names(o);
    else
        names = sort(names);
    end
end

function g = read_gray(pathname)
% Read one frame as grayscale double in [0..1] whatever the file's bit
% depth: integer classes are scaled by their OWN intmax (uint8 -> /255,
% uint16 -> /65535, ...), so 10/12/16-bit scientific-camera series satisfy
% the engine's [0..1] contract (P.black_thr, quality gates). A fixed /255
% would leave 16-bit data ~257x out of range. Float images already above 1
% are scaled by their max. Same as /255 for 8-bit input.
    I = imread(pathname);
    if isinteger(I), s = double(intmax(class(I))); else, s = 1; end
    g = to_gray(I) / s;
    mx = max(g(:));
    if mx > 1, g = g / mx; end
end

function roi = detect_band(img, W, H)
% ROI (bounding box) around the specimen band. It feeds BOTH the coarse
% large-motion frame-to-frame shift search AND the grid extent.
%   ROWS  - the tight textured band (texture > mean). On a LARGE-motion series
%           (hundreds of px per step) the coarse search only locks onto the
%           real displacement when the tracked strip is the high-contrast
%           specimen CORE; broadening the rows pulls in low-motion background
%           and collapses the first step ("no step completed").
%   COLS  - FULL WIDTH, always. Thresholding columns cropped a wide beam to its
%           high-contrast half ("only half the specimen measured"); the beam
%           spans the frame, so columns are never trimmed.
% For a specimen that does NOT span the width (or any odd framing) draw a
% manual ROI (Select ROI) - it bypasses this and is used exactly as drawn.
% Decimated for speed; returns [x y w h] px.
    dd = 4;
    A = img(1:dd:end, 1:dd:end);
    rs = std(A, 0, 2);                 % per-row texture (across columns)
    rows = find(rs > mean(rs));
    if isempty(rows), roi = [1 1 W-1 H-1]; return; end
    pad = 3;
    r_lo = max(1, (min(rows) - pad) * dd);
    r_hi = min(H, (max(rows) + pad) * dd);
    roi = [1, r_lo, W - 1, r_hi - r_lo];   % full width, tracked band in rows
end

function xy = best_texture_xy(img, roi, W, H)
% Point of maximum local texture (std) inside the ROI - the safest place
% to put the DIC seed. Decimated, box statistics via conv2.
    dd = 4; w = 11;
    A = img(1:dd:end, 1:dd:end);
    K = ones(w) / w^2;
    m = conv2(A, K, 'same');
    v = conv2(A.^2, K, 'same') - m.^2;
    [Ha, Wa] = size(A);
    r1 = min(max(1, ceil(roi(2)/dd) + w), Ha);
    r2 = min(max(1, floor((roi(2)+roi(4))/dd) - w), Ha);
    c1 = min(max(1, ceil(roi(1)/dd) + w), Wa);
    c2 = min(max(1, floor((roi(1)+roi(3))/dd) - w), Wa);
    if r2 <= r1 || c2 <= c1, xy = []; return; end
    sub = v(r1:r2, c1:c2);
    [~, im] = max(sub(:));
    [ri, ci] = ind2sub(size(sub), im);
    xy = [min(W, (c1 + ci - 1) * dd), min(H, (r1 + ri - 1) * dd)];
end

function [sx, sy, best] = coarse_shift(prev, cur, roi, ds, sRr, sRc, prior)
% Integer shift of the specimen between consecutive frames: decimated
% ZNCC match of a template taken from the ROI center of `prev`. The
% search window is CENTERED on the previous step's shift (temporal prior).
    A = prev(1:ds:end, 1:ds:end);
    B = cur (1:ds:end, 1:ds:end);
    [Ha, Wa] = size(A);
    r_lo = max(1, floor(roi(2)/ds));  r_hi = min(Ha, ceil((roi(2)+roi(4))/ds));
    tr = round((r_lo + r_hi) / 2);
    th = max(8, round((r_hi - r_lo) * 0.4));
    tw = max(8, round(Wa * 0.25));
    r0 = min(max(1, tr - round(th/2)), max(1, Ha - th));
    c0 = min(max(1, round(Wa/2 - tw/2)), max(1, Wa - tw));
    T = A(r0:r0+th-1, c0:c0+tw-1);
    Tz = T - mean(T(:));  Tn = norm(Tz(:));
    p_r = round(prior(2) / ds);  p_c = round(prior(1) / ds);
    best = -2; dr = p_r; dc = p_c;
    if Tn < eps, sx = prior(1); sy = prior(2); return; end
    for ddr = p_r-sRr : p_r+sRr
        r = r0 + ddr;
        if r < 1 || r + th - 1 > Ha, continue; end
        for ddc = p_c-sRc : p_c+sRc
            c = c0 + ddc;
            if c < 1 || c + tw - 1 > Wa, continue; end
            Pp = B(r:r+th-1, c:c+tw-1);
            Pz = Pp - mean(Pp(:));  Pn = norm(Pz(:));
            if Pn < eps, continue; end
            v = sum(sum(Tz .* Pz)) / (Tn * Pn);
            if v > best, best = v; dr = ddr; dc = ddc; end
        end
    end
    sx = dc * ds;  sy = dr * ds;
end

function dump_pair(prev, cur, roi, s, outDir, k)
% Save decimated ROI crops of the failing pair for eyeball diagnosis:
% prev ROI, and the same window in cur shifted by the coarse estimate.
    [H, W] = size(prev);
    r1 = max(1, roi(2));  r2 = min(H, roi(2) + roi(4));
    c1 = max(1, roi(1));  c2 = min(W, roi(1) + roi(3));
    Ap = prev(r1:4:r2, c1:4:c2);
    rs = max(1, r1 + s(2));  re = min(H, r2 + s(2));
    cs = max(1, c1 + s(1));  ce = min(W, c2 + s(1));
    Ac = cur(rs:4:re, cs:4:ce);
    imwrite(norm01(Ap), fullfile(outDir, sprintf('debug_step%02d_prev.png', k)));
    imwrite(norm01(Ac), fullfile(outDir, sprintf('debug_step%02d_cur.png', k)));
end

function g = norm01(A)
    lo = min(A(:)); hi = max(A(:));
    if hi > lo, g = (A - lo) / (hi - lo); else, g = zeros(size(A)); end
end

function r = clamp_roi(roi, W, H)
    x = min(max(1, roi(1)), W - 2);
    y = min(max(1, roi(2)), H - 2);
    r = [x, y, min(roi(3), W - x - 1), min(roi(4), H - y - 1)];
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end

function g = to_gray(img)
    if ndims(img) == 3 && size(img, 3) >= 3
        img = double(img);
        g = 0.2989*img(:,:,1) + 0.5870*img(:,:,2) + 0.1140*img(:,:,3);
    else
        g = double(img);
    end
end

function K = gk(sigma)
    r = ceil(3 * sigma);
    [x, y] = meshgrid(-r:r);
    K = exp(-(x.^2 + y.^2) / (2 * sigma^2));
    K = K / sum(K(:));
end
