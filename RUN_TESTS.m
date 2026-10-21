function RUN_TESTS(mode, dataDir)
% RUN_TESTS - automated regression harness: runs every analysis and checks
% the key numbers against reference intervals.
%
% Usage from the project folder:
%       RUN_TESTS                  % main set, ~1-2 min
%       RUN_TESTS('stereo-data')   % + welded-plate stereo vs Vic-3D
%                                  %   (OPT-IN: needs the external dataset,
%                                  %    takes minutes; base suite is
%                                  %    independent of any external data)
%       RUN_TESTS('stereo-data', '<path to DIC 3D dataset>')  % custom path
%
% The report is printed and saved to TEST_REPORT.txt.

if nargin < 1, mode = 'base'; end
if nargin < 2 || isempty(dataDir)
    % Portable default for the opt-in stereo-data set (~2.5 GB, not in git):
    % (1) the LITEDIC_STEREO_DATA environment variable, else (2) a few common
    % local locations. No user-specific absolute path is baked into source;
    % if none resolves, test_stereo_dataset prints how to pass the path.
    rootDir = fileparts(mfilename('fullpath'));
    if ispc, homeDir = getenv('USERPROFILE'); else, homeDir = getenv('HOME'); end
    dataDir = getenv('LITEDIC_STEREO_DATA');
    if isempty(dataDir)
        cand = { fullfile(homeDir, 'DIC', 'DIC 3D dataset'), ...
                 fullfile(fileparts(rootDir), 'LiteDIC app', 'examples', 'DIC 3D dataset'), ...
                 fullfile(rootDir, '12_Datasets', 'DIC 3D dataset') };
        for ii = 1:numel(cand)
            if isfolder(cand{ii}), dataDir = cand{ii}; break; end
        end
        if isempty(dataDir), dataDir = cand{1}; end   % name a sensible path in the error
    end
end
do_stereo_data = strcmpi(mode, 'stereo-data');
root = fileparts(mfilename('fullpath'));
t_all = tic;
G = {};   % report rows: {group, check name, ok(0/1/-1), details}

% ---------- preflight ----------
% rehash: refresh MATLAB's path/file cache - REQUIRED after any file
% operations done outside MATLAB (sync tools, file managers, OS commands),
% otherwise run() may use a stale cached classification of the files.
rehash;
ensure_images(root);   % verify the DIC image series is present
file_check(root);      % print sizes + first lines as MATLAB sees the files

G = run_case(G, 'DIC coarse (dic_coarse.m)',    @test_corr);
G = run_case(G, 'DIC subpixel (dic_icgn.m)',    @test_dicpro);
G = run_case(G, 'DIC batch kernel (dic_engine_batch.m)', @test_batch);
G = run_case(G, 'DIC bias S-curve (bias_scurve.m)', @test_scurve);
G = run_case(G, 'DIC PLS strain vs grid step (pls_strain.m)', @test_grid_refine);
G = run_case(G, 'DIC large-shift seed (run_series_2d path)', @test_series2d_bigshift);
G = run_case(G, 'DIC real-data inputs (bit depth, frame order, size guard)', @test_realdata_inputs);
G = run_case(G, 'DIC engine AOI mask hook (dic_engine P.mask)', @test_engine_mask);
G = run_case(G, 'DIC multi-seed barrier crossing (dic_engine)', @test_multiseed);
G = run_case(G, 'DIC per-point uncertainty (dic_engine R.SU/R.SV)', @test_sigma);
G = run_case(G, 'DIC virtual extensometer (virtual_extensometer.m)', @test_extensometer);
G = run_case(G, 'DIC mm/px scale (dic_scale.m)', @test_scale);
G = run_case(G, 'DIC vs DICe gold (validate_dice.m)', @test_dice_gold);
G = run_case(G, 'Stereo intrinsics (11_Stereo)', @test_stereo_cal);
G = run_case(G, 'Stereo extrinsics+triangulation (11_Stereo)', @test_stereo3);
G = run_case(G, 'Stereo DIC vs DICe gold (11_Stereo)', @test_stereo_gold);
G = run_case(G, 'Surface strain (surf_strain.m)', @test_surf_strain);
G = run_case(G, 'GT dataset factory (gen_dataset.m)', @test_gt_gen);
G = run_case(G, 'External validator (13_Validator)', @test_validator);
if do_stereo_data
    G = run_case(G, 'Stereo dataset vs Vic-3D (frozen, opt-in)', ...
        @() test_stereo_dataset(dataDir));
end

% ---------- report ----------
npass = sum(cellfun(@(x) isequal(x,1),  G(:,3)));
nfail = sum(cellfun(@(x) isequal(x,0),  G(:,3)));
nerr  = sum(cellfun(@(x) isequal(x,-1), G(:,3)));
L = {};
L{end+1} = sprintf('PROJECT TEST REPORT - %s', char(datetime('now','Format','yyyy-MM-dd HH:mm')));
L{end+1} = sprintf('MATLAB %s | mode: %s | time: %.0f s', version, mode, toc(t_all));
L{end+1} = sprintf('TOTAL: %d PASS, %d FAIL, %d ERRORS', npass, nfail, nerr);
L{end+1} = '';
grp = '';
for k = 1:size(G,1)
    if ~strcmp(G{k,1}, grp), grp = G{k,1}; L{end+1} = ['--- ' grp]; end %#ok<AGROW>
    switch G{k,3}
        case 1,  tag = '[PASS]';
        case 0,  tag = '[FAIL]';
        otherwise, tag = '[ERROR]';
    end
    L{end+1} = sprintf('%s %s%s', tag, G{k,2}, G{k,4}); %#ok<AGROW>
end
rep = strjoin(L, newline);
fprintf('\n%s\n', rep);
fid = fopen(fullfile(root, 'TEST_REPORT.txt'), 'w', 'n', 'UTF-8');
if fid > 0, fwrite(fid, rep); fclose(fid); end
fprintf('\nReport saved to TEST_REPORT.txt.\n');
end

% ===================== framework =====================

function G = run_case(G, name, fcn)
    fprintf('>> %-46s ', [name ' ...']);
    t0 = tic;
    try
        rows = fcn();
        for k = 1:size(rows,1)
            G(end+1,:) = {name, rows{k,1}, rows{k,2}, rows{k,3}}; %#ok<AGROW>
        end
        nf = sum(cellfun(@(x) isequal(x,0), rows(:,2)));
        if nf == 0, fprintf('OK (%.0f s)\n', toc(t0));
        else,       fprintf('%d FAIL (%.0f s)\n', nf, toc(t0));
        end
    catch err
        G(end+1,:) = {name, 'execution', -1, [': ' err.message]}; %#ok<AGROW>
        fprintf('ERROR: %s\n', err.message);
    end
    close all;
end

function rows = chk(rows, name, val, lo, hi, fmt)
% Check that val falls inside the reference interval [lo, hi].
    if nargin < 6, fmt = '%.4g'; end
    ok = isnumeric(val) && isscalar(val) && isfinite(val) && val >= lo && val <= hi;
    rows(end+1,:) = {name, double(ok), ...
        sprintf([': ' fmt ' (reference %.4g .. %.4g)'], val, lo, hi)};
end

function p = proj(varargin)
    p = fullfile(fileparts(mfilename('fullpath')), varargin{:});
end

function ensure_images(root)
% Verify the DIC image series is present and non-empty. The series is
% discovered by pattern (same convention as the DIC scripts themselves).
    dst = fullfile(root, '6_DIC');
    d_ref = dir(fullfile(dst, 'plate_ref.tiff'));
    d_def = dir(fullfile(dst, 'plate_def_*.tiff'));
    d_def = d_def([d_def.bytes] > 0);
    if isempty(d_ref) || d_ref(1).bytes == 0
        warning('DIC reference frame plate_ref.tiff is missing or empty in %s.', dst);
    end
    if isempty(d_def)
        warning('No DIC deformed frames (plate_def_*.tiff) found in %s.', dst);
    else
        fprintf('DIC series: plate_ref.tiff + %d deformed frames.\n', numel(d_def));
    end
end

function file_check(root)
% Diagnostic listing: size and first line of every script AS MATLAB READS
% IT from disk. If a file shows 0 bytes or garbage here, the file on disk
% is broken even if other tools claim it is fine.
    files = { ...
        fullfile('6_DIC','dic_coarse.m'), ...
        fullfile('6_DIC','dic_icgn.m'), ...
        fullfile('6_DIC','dic_engine.m'), ...
        fullfile('6_DIC','dic_engine_batch.m'), ...
        fullfile('6_DIC','pls_strain.m'), ...
        fullfile('6_DIC','virtual_extensometer.m'), ...
        fullfile('6_DIC','dic_scale.m'), ...
        fullfile('6_DIC','scale_from_length.m'), ...
        fullfile('6_DIC','plate_ref.tiff'), ...
        fullfile('9_Validation_DICe','validate_dice.m'), ...
        fullfile('9_Validation_DICe','dice_reference','DICe_solution_0.txt'), ...
        fullfile('10_Challenge2','challenge2_benchmark.m'), ...
        'START.m'};
    fprintf('--- File check (as MATLAB sees the disk) ---\n');
    for k = 1:numel(files)
        p = fullfile(root, files{k});
        d = dir(p);
        if isempty(d)
            fprintf('%-45s MISSING\n', files{k});
            continue;
        end
        info = sprintf('%7d B', d(1).bytes);
        if endsWith(p, '.m')
            fid = fopen(p, 'r');
            l1 = fgetl(fid);
            fclose(fid);
            if ~ischar(l1), l1 = '<empty>'; end
            info = [info '  | ' l1(1:min(end,48))]; %#ok<AGROW>
        end
        fprintf('%-45s %s\n', files{k}, info);
    end
    fprintf('--------------------------------------------\n\n');
end

function v = prc(x, q)
% Percentile without the Statistics Toolbox.
    x = sort(x(isfinite(x))); n = numel(x);
    v = x(min(n, max(1, round(q/100*n))));
end

% ===================== tests =====================

function rows = test_corr()
% dic_coarse runs on plate_ref -> last frame of the series.
% Regression intervals calibrated on the canonical TIFF series
% (plate_ref + def_01..11), R2025b, 2026-06-12: dy -8.51, |dx| 0.68,
% NCC 0.960. Recalibrate if the image series changes.
    evalc('run(proj(''6_DIC'', ''dic_coarse.m''))');
    rows = {};
    rows = chk(rows, 'mean dy (downward load), px', mean(displacement_y(:), 'omitnan'), -9.0, -8.0, '%.2f');
    rows = chk(rows, 'mean |dx| << |dy|, px', abs(mean(displacement_x(:), 'omitnan')), 0.4, 1.0, '%.2f');
    rows = chk(rows, 'mean NCC', mean(ncc_values(:), 'omitnan'), 0.945, 0.975, '%.3f');
end

function rows = test_dicpro()
    evalc('run(proj(''6_DIC'', ''dic_icgn.m''))');
    rows = {};
    % self-test (variables eu, ev, smax persist after the script)
    rows = chk(rows, 'self-test: |error u|, px', abs(eu), 0, 0.02, '%.5f');
    rows = chk(rows, 'self-test: |error v|, px', abs(ev), 0, 0.02, '%.5f');
    rows = chk(rows, 'self-test: spurious strains', smax, 0, 5e-4, '%.1e');
    % Loading series ("series" shadows a toolbox function - copy first).
    % These are the README checkpoints as INVARIANTS of any monotonic
    % loading series, not values pinned to specific frames.
    % Regression intervals calibrated on the canonical TIFF series
    % (plate_ref + def_01..11), R2025b, 2026-06-12: last |v| 8.50,
    % mean ZNCC 0.9976, max|eyy| 0.056, max/p99 2.30. Recalibrate if the
    % image series changes; the monotonicity check stays series-agnostic.
    sv = [series.mean_v]; sz = [series.zncc];
    rows = chk(rows, 'series length, frames', numel(sv), 2, Inf, '%.0f');
    rows = chk(rows, '|v| grows monotonically frame to frame', ...
        double(all(diff(abs(sv)) > 0)), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'last frame: |v|, px', abs(sv(end)), 8.3, 8.7, '%.2f');
    rows = chk(rows, 'mean ZNCC of the series', mean(sz), 0.995, 0.999, '%.4f');
    % strains of the last frame: concentration at the hole edge, no fake
    % peaks (max should be a few times the p99 level, not orders above)
    e_max = max(abs(eyy(:)), [], 'omitnan');
    e_p99 = prc(abs(eyy), 99);
    rows = chk(rows, 'max |eyy| (edge concentration)', e_max, 0.045, 0.070, '%.3f');
    rows = chk(rows, 'concentration ratio max/p99 |eyy|', e_max / e_p99, 1.9, 2.8, '%.2f');
end

function rows = test_grid_refine()
% Grid-step behaviour of the PLS strain field. The strain window is given
% in PIXELS, so refining the analysis grid is supposed to leave the
% physical smoothing length - and hence the strain field - unchanged.
% That holds only up to the rounding of the half-window to a whole number
% of nodes (pls_strain: h = max(1, round(win_px/(2*step)))), and both
% halves of the statement are checked here: exact invariance where the
% rounding is exact, and a bounded departure where it is not.
%
% The field is analytic, so no correlation is involved and the check is
% about pls_strain alone: v = A*sin(2*pi*y/lam) has
% eyy = (2*pi*A/lam)*cos(2*pi*y/lam), and a plane fit over a window w
% attenuates that amplitude by T(w) = 3*(sin q - q*cos q)/q^3 with
% q = pi*w/lam - a function of w in PIXELS only.
    addpath(proj('6_DIC'));
    lam = 120; A = 0.5; L = 600;                 % px
    amp0 = 2*pi*A/lam;                           % true eyy amplitude
    rows = {};

    % (a) the implementation reproduces the exact transfer of its own
    %     estimator, for a fine grid, a coarse one, and a rounded window
    [y1, e1] = sine_eyy(4,  40, lam, A, L);      % h = 5, span 40 px
    [y2, e2] = sine_eyy(10, 40, lam, A, L);      % h = 2, span 40 px
    [y3, e3] = sine_eyy(8,  40, lam, A, L);      % h = 3, span 48 px (rounded)
    a1 = harm_amp(y1, e1, lam);
    a2 = harm_amp(y2, e2, lam);
    a3 = harm_amp(y3, e3, lam);
    rows = chk(rows, 'step 4  (h=5): measured / exact transfer', ...
        a1 / (amp0 * pls_transfer(4, 5, lam)),  0.995, 1.005, '%.4f');
    rows = chk(rows, 'step 10 (h=2): measured / exact transfer', ...
        a2 / (amp0 * pls_transfer(10, 2, lam)), 0.995, 1.005, '%.4f');
    rows = chk(rows, 'step 8  (h=3): measured / exact transfer', ...
        a3 / (amp0 * pls_transfer(8, 3, lam)),  0.995, 1.005, '%.4f');

    % (b) invariance with the window fixed in PIXELS. Steps 4 and 10 give
    %     the same 40 px span, yet the fields still differ by ~2.7 %: the
    %     plane fit is a discrete sum over the nodes inside the window, so
    %     a coarser grid changes the quadrature even at a fixed span. This
    %     is where the invariance breaks, see pls_strain.m.
    [dmax, ncom] = eyy_common_diff(y1, e1, y2, e2);
    rows = chk(rows, 'window in px: common nodes (step 4 vs 10)', ...
        ncom, 20, Inf, '%.0f');
    rows = chk(rows, 'window in px: departure, step 4 vs 10', ...
        abs(a1/a2 - 1), 0.02, 0.035, '%.4f');
    rows = chk(rows, 'window in px: max |d eyy| / amplitude', ...
        dmax / amp0, 0.015, 0.035, '%.4f');

    % (c) the same comparison with the window fixed in NODES instead
    %     (h = 2 on both grids, i.e. spans of 16 and 40 px): about five
    %     times larger departure.
    [~, e4] = sine_eyy(4, 16, lam, A, L);        % h = 2, span 16 px
    a4 = harm_amp(y1, e4, lam);
    rows = chk(rows, 'window in nodes: departure, step 4 vs 10', ...
        abs(a4/a2 - 1), 0.12, 0.17, '%.4f');
end

function [gy, eyy] = sine_eyy(step, win_px, lam, A, L)
% eyy field of v = A*sin(2*pi*y/lam) on a grid of the given step, keeping
% only the interior band where every strain window is complete.
    gy = (0:step:L)';
    V  = A * sin(2*pi*gy/lam) * ones(1, numel(0:step:L));
    [~, eyy] = pls_strain(zeros(size(V)), V, step, win_px, false);
    eyy = eyy(:, ceil(size(eyy,2)/2));           % field varies in y only
    keep = gy >= 80 & gy <= L - 80;
    gy = gy(keep); eyy = eyy(keep);
end

function [dmax, n] = eyy_common_diff(y1, e1, y2, e2)
% Largest |difference| between two strain fields at the nodes they share.
    [~, i1, i2] = intersect(round(y1, 6), round(y2, 6));
    d = e1(i1) - e2(i2);
    d = d(isfinite(d));
    n = numel(d);
    if n == 0, dmax = Inf; else, dmax = max(abs(d)); end
end

function a = harm_amp(y, e, lam)
% First-harmonic amplitude of e(y) at the known period - grid-independent,
% unlike max(|e|), which depends on where the nodes fall on the cosine.
    m = isfinite(e); yy = y(m); ee = e(m);
    c = [cos(2*pi*yy/lam), sin(2*pi*yy/lam), ones(size(yy))] \ ee;
    a = hypot(c(1), c(2));
end

function t = pls_transfer(step, h, lam)
% Exact amplitude transfer of the PLS slope estimator on a grid: the fit
% is the discrete sum slope = sum(x.*q)/sum(x.^2) over the 2h+1 nodes of
% the window, so the transfer depends on the node spacing as well as on
% the span 2*h*step. (The continuous limit is 3*(sin q - q cos q)/q^3.)
    x = (-h:h)' * step;
    t = sum(x .* sin(2*pi*x/lam)) / ((2*pi/lam) * sum(x.^2));
end

function rows = test_batch()
% The batched IC-GN kernel must (a) recover a known rigid subpixel shift
% and (b) agree with the sequential dic_engine to within the convergence
% tolerance scale (both converge to the same optimum, conv_tol = 1e-3).
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(220), fspecial_gauss(2.0), 'same');
    img = img(11:210, 11:210);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 0.30; TV = -0.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    P = struct('subset_size', 41, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3);
    Rb = dic_engine_batch(img, shifted, P);
    Rs = dic_engine(img, shifted, P);
    m = isfinite(Rb.U) & isfinite(Rs.U);
    dmax = max(max(abs(Rb.U(m) - Rs.U(m))), max(abs(Rb.V(m) - Rs.V(m))));
    rows = {};
    rows = chk(rows, 'batch: |error u|, px', abs(mean(Rb.U(:), 'omitnan') - TU), 0, 0.02, '%.5f');
    rows = chk(rows, 'batch: |error v|, px', abs(mean(Rb.V(:), 'omitnan') - TV), 0, 0.02, '%.5f');
    rows = chk(rows, 'batch vs sequential: common points', nnz(m), 100, Inf, '%.0f');
    rows = chk(rows, 'batch vs sequential: max |diff|, px', dmax, 0, 0.005, '%.2e');
end

function rows = test_scurve()
% Mini version of the 6_DIC/bias_scurve.m sweep: a periodic speckle is
% Fourier-phase-shifted (generator adds no interpolation bias of its own)
% by a few subpixel amounts; the engine's mean error must stay at the
% interpolation-bias scale and there must be no cross-talk into v.
% The full 21-shift, 2-pattern, 2-subset S-curve lives in bias_scurve.m.
    addpath(proj('6_DIC'));
    rng(0);
    N = 128;
    W = fft2(rand(N));
    k = ifftshift(-floor(N/2):ceil(N/2)-1) * (2*pi/N);
    [KX, KY] = meshgrid(k, k);
    img = real(ifft2(W .* exp(-(KX.^2 + KY.^2) * 2.0^2 / 2)));   % sigma 2.0
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    P = struct('subset_size', 21, 'grid_step', 4, 'search_range', 2, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-4);
    deltas = [0.10 0.25 0.40];
    bu = nan(size(deltas)); bv = bu; np = bu;
    for i = 1:numel(deltas)
        shifted = real(ifft(fft(img, [], 2) .* exp(-1i * deltas(i) * k), [], 2));
        R = dic_engine_batch(img, shifted, P);
        u = R.U(isfinite(R.U));  v = R.V(isfinite(R.V));
        bu(i) = mean(u) - deltas(i);
        bv(i) = mean(v);
        np(i) = numel(u);
    end
    rows = {};
    rows = chk(rows, 'S-curve: max |u bias|, px', max(abs(bu)), 0, 0.02, '%.2e');
    rows = chk(rows, 'S-curve: max |v cross-talk|, px', max(abs(bv)), 0, 0.005, '%.2e');
    rows = chk(rows, 'S-curve: points per shift', min(np), 300, Inf, '%.0f');
end

function rows = test_scale()
% Optional physical scale (mm/px) is a pure post-processing layer on the
% engine output (dic_scale.m). It must (a) reproduce the pixel fields
% bit-for-bit at unit scale, (b) multiply displacements AND grid
% coordinates by mm_per_px at any other scale, and (c) leave the
% dimensionless strain field untouched. scale_from_length.m recovers
% mm_per_px from a known real-world length between two points. Every check
% here is an exact logical invariant (bit-exact IEEE multiply / deadband
% comparison), so no image-dependent intervals need pinning.
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(220), fspecial_gauss(2.0), 'same');
    img = img(11:210, 11:210);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 0.30; TV = -0.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    P = struct('subset_size', 41, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3);
    R = dic_engine_batch(img, shifted, P);
    [R.exx, R.eyy, R.exy] = pls_strain(R.U, R.V, P.grid_step, 30, false);

    % (a) unit scale: every numeric field identical (NaNs included)
    R1 = dic_scale(R, 1.0);
    unit_ok = isequaln(R1.U, R.U) && isequaln(R1.V, R.V) && ...
              isequaln(R1.gx, R.gx) && isequaln(R1.gy, R.gy) && ...
              isequaln(R1.eyy, R.eyy);

    % (b) scale 0.05: displacements and coordinates x0.05, exactly
    s = 0.05;
    Rs = dic_scale(R, s);
    disp_ok  = isequaln(Rs.U, s*R.U) && isequaln(Rs.V, s*R.V);
    coord_ok = isequaln(Rs.gx, s*R.gx) && isequaln(Rs.gy, s*R.gy);

    % (c) strains are dimensionless: unchanged at any scale
    strain_ok = isequaln(Rs.exx, R.exx) && isequaln(Rs.eyy, R.eyy) && ...
                isequaln(Rs.exy, R.exy);

    % helper: 25.4 mm across a 512 px span -> exact mm/px
    mmpp = scale_from_length(25.4, [100 300], [612 300]);
    helper_ok = abs(mmpp - 25.4/512) < 1e-12;

    rows = {};
    rows = chk(rows, 'scale=1: length + strain fields bit-for-bit', double(unit_ok), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'scale=0.05: U,V multiplied exactly', double(disp_ok), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'scale=0.05: gx,gy multiplied exactly', double(coord_ok), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'strain field unchanged under scaling', double(strain_ok), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'units tag set to mm', double(strcmp(Rs.units, 'mm')), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'scale_from_length: 25.4 mm / 512 px correct', double(helper_ok), 0.5, 1.5, '%.0f');
end

function K = fspecial_gauss(sigma)
% Gaussian kernel (same recipe as the dic_icgn.m self-test; no toolboxes).
    r = ceil(3*sigma);
    [x, y] = meshgrid(-r:r);
    K = exp(-(x.^2 + y.^2) / (2*sigma^2));
    K = K / sum(K(:));
end

function rows = test_series2d_bigshift()
% Large rigid-shift recovery through P.seed_shift - needed by the
% run_series_2d front end. The +/-search_range
% (15 px) integer seed search is CENTERED on a coarse integer prior, so a
% shift of tens-to-hundreds of px per step is recovered without enlarging
% the search window. This is exactly the synthetic self-test embedded in
% run_series_2d.m (rigid shift +60.30/-40.20 px, coarse prior [60 -40]),
% wired in here as an always-on regression so the seed_shift path cannot
% silently break. Deterministic (rng 0), toolbox-free, no external data.
% Bounds mirror the driver's own accept threshold (|err| < 0.05 px); the
% 400x400 frame at grid step 10 (margin 37) yields a 33x33 node grid.
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(420), fspecial_gauss(2.0), 'same');
    img = img(11:410, 11:410);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 60.30; TV = -40.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    P = struct('subset_size', 41, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3, ...
        'seed_shift', [60 -40]);
    R = dic_engine(img, shifted, P);
    eu = mean(R.U(:), 'omitnan') - TU;
    ev = mean(R.V(:), 'omitnan') - TV;
    rows = {};
    rows = chk(rows, 'seed_shift: |error u|, px', abs(eu), 0, 0.05, '%.5f');
    rows = chk(rows, 'seed_shift: |error v|, px', abs(ev), 0, 0.05, '%.5f');
    rows = chk(rows, 'seed_shift: valid points (no collapse)', nnz(~isnan(R.U)), 500, Inf, '%.0f');
end

function rows = test_realdata_inputs()
% Real-data input checks for the 2D series driver:
% (1) non-zero-padded frame numbers are processed in natural
%     (chronological) order - a lexicographic sort reads frame_10 before
%     frame_2 and silently scrambles the series;
% (2) frames are normalized by their OWN bit depth: a 16-bit series must
%     give the same fields AND the same black-region (dark hole) mask as
%     its 8-bit twin (a hardcoded /255 puts 16-bit data ~257x out of
%     range and the hole mask stops working);
% (3) a stray mixed-size file dies with a clear message, not a cryptic
%     interp error;
% (4) the driver's built-in self-test is PINNED to its own subset/step: a
%     fine-speckle user subset (15 px) must not fail it spuriously.
% All synthetic, deterministic (rng 0), toolbox-free, no external data.
    addpath(proj('6_DIC'));
    rng(0);
    M = conv2(rand(400, 440), fspecial_gauss(2.0), 'same');
    M = M(11:390, 11:430);                        % 380 x 420 master scene
    M = (M - min(M(:))) / (max(M(:)) - min(M(:)));
    M = 0.05 + 0.9 * M;                           % keep off the int limits
    [xx, yy] = meshgrid(1:size(M,2), 1:size(M,1));
    M(hypot(xx - 150, yy - 170) < 35) = 0.02;     % dark hole (black-region mask)
    dy = 3;                                       % rigid shift per frame, px
    crop = @(k) M((1:300) + (k-1)*dy, 1:340);     % scene moves up: v = -dy/step
    args = {'SelfTest', false, 'Frames', false, 'MaxShift', [8 8], ...
            'GridStep', 10, 'SubsetSize', 41};
    rows = {};

    % ---- (1) natural frame order + the 8-bit reference run ----
    d8 = [tempname 'a']; mkdir(d8);
    cu8 = onCleanup(@() rmdir(d8, 's')); %#ok<NASGU>
    nm = {'frame_1.png', 'frame_2.png', 'frame_10.png'};   % 10 sorts after 2
    for k = 1:3
        imwrite(uint8(round(crop(k) * 255)), fullfile(d8, nm{k}));
    end
    evalc('run_series_2d(d8, args{:})');
    S8 = load(fullfile(d8, 'litedic_out', 'series_2d_results.mat'), ...
        'names', 'cumV', 'alive', 'X0', 'Y0');
    rows = chk(rows, 'non-padded frame names in natural order', ...
        double(isequal(S8.names(:).', nm)), 0.5, 1.5, '%.0f');
    cv8 = mean(S8.cumV(S8.alive), 'omitnan');
    rows = chk(rows, '8-bit run: cum v after 2 steps of -3 px, px', ...
        cv8, -2*dy - 0.05, -2*dy + 0.05, '%.3f');

    % ---- (2) 16-bit twin: same fields, same hole mask ----
    d16 = [tempname 'b']; mkdir(d16);
    cu16 = onCleanup(@() rmdir(d16, 's')); %#ok<NASGU>
    for k = 1:3
        imwrite(uint16(round(crop(k) * 65535)), fullfile(d16, sprintf('frame_%02d.png', k)));
    end
    evalc('run_series_2d(d16, args{:})');
    S16 = load(fullfile(d16, 'litedic_out', 'series_2d_results.mat'), ...
        'cumV', 'alive');
    inHole = hypot(S8.X0 - 150, S8.Y0 - 170) < 8;  % nodes deep inside the hole
    nH8  = nnz(inHole & ~isfinite(S8.cumV));
    nH16 = nnz(inHole & ~isfinite(S16.cumV));
    rows = chk(rows, 'hole nodes masked in the 8-bit run', nH8, 1, Inf, '%.0f');
    rows = chk(rows, '16-bit run: hole nodes masked like 8-bit', ...
        nH16, nH8 - 0.5, nH8 + 0.5, '%.0f');
    rows = chk(rows, '16-bit vs 8-bit: differing valid nodes', ...
        nnz(xor(isfinite(S8.cumV), isfinite(S16.cumV))), 0, 2, '%.0f');
    both = S8.alive & S16.alive;
    rows = chk(rows, '16-bit vs 8-bit: max |cum v diff|, px', ...
        max(abs(S8.cumV(both) - S16.cumV(both))), 0, 0.02, '%.4f');

    % ---- (3) mixed-size stray file -> clear error ----
    dm = [tempname 'c']; mkdir(dm);
    cum = onCleanup(@() rmdir(dm, 's')); %#ok<NASGU>
    imwrite(uint8(round(crop(1) * 255)), fullfile(dm, 'frame_01.png'));
    imwrite(uint8(round(crop(2) * 255)), fullfile(dm, 'frame_02.png'));
    C3 = crop(3);
    imwrite(uint8(round(C3(1:150, 1:170) * 255)), fullfile(dm, 'frame_03.png'));
    okmsg = false;
    try
        evalc('run_series_2d(dm, args{:})');
    catch err
        okmsg = contains(err.message, 'mixed-size');
    end
    rows = chk(rows, 'mixed-size stray frame errors clearly', ...
        double(okmsg), 0.5, 1.5, '%.0f');

    % ---- (4) pinned self-test: fine-speckle subset must not fail it ----
    ok15 = true;
    try
        evalc(['run_series_2d(d16, ''SelfTest'', true, ''Frames'', false, ' ...
               '''MaxShift'', [8 8], ''GridStep'', 12, ''SubsetSize'', 15)']);
    catch
        ok15 = false;
    end
    rows = chk(rows, 'pinned self-test: SubsetSize 15 run survives', ...
        double(ok15), 0.5, 1.5, '%.0f');
end

function rows = test_engine_mask()
% Optional AOI mask hook (dic_engine P.mask) - the engine half of the
% polygon-AOI feature. A full-frame boolean restricts the grid alongside
% P.roi so fixtures/glare can be excluded. Every check is an exact logical
% invariant (no image-dependent interval to pin), deterministic (rng 0),
% toolbox-free, no external data:
%   (a) an all-true mask reproduces the no-mask fields BIT-FOR-BIT (the mask
%       path must be inert when it excludes nothing);
%   (b) nodes in a masked-out band are all NaN (the mask actually excludes);
%   (c) the kept nodes match the no-mask run to convergence tolerance (only
%       the boundary warm-start can differ; both reach the same optimum).
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(220), fspecial_gauss(2.0), 'same');
    img = img(11:210, 11:210);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 0.30; TV = -0.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    P = struct('subset_size', 41, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3);
    R0 = dic_engine(img, shifted, P);

    % (a) all-true mask must be identical to no mask
    Pm = P; Pm.mask = true(size(img));
    R1 = dic_engine(img, shifted, Pm);
    same = isequaln(R1.U, R0.U) && isequaln(R1.V, R0.V);

    % (b)+(c) exclude the bottom band of the frame (image rows > 110)
    m = true(size(img)); m(111:end, :) = false;
    Pm2 = P; Pm2.mask = m;
    R2 = dic_engine(img, shifted, Pm2);
    inband  = R2.gy(:) <= 110;          % grid rows kept by the mask
    outband = R2.gy(:) >  110;          % grid rows excluded by the mask
    Uout = R2.U(outband, :);
    allnan_out = ~any(isfinite(Uout(:)));
    U2 = R2.U(inband, :); U0 = R0.U(inband, :);
    V2 = R2.V(inband, :); V0 = R0.V(inband, :);
    mm = isfinite(U2) & isfinite(U0);
    dmax = max(max(abs(U2(mm) - U0(mm))), max(abs(V2(mm) - V0(mm))));

    rows = {};
    rows = chk(rows, 'all-true mask == no mask (bit-for-bit)', double(same), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'masked-out nodes are all NaN', double(allnan_out), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'kept nodes match no-mask run, max|diff| px', dmax, 0, 0.005, '%.2e');
    rows = chk(rows, 'kept nodes stay valid', nnz(mm), 30, Inf, '%.0f');
end

function rows = test_sigma()
% Per-point uncertainty (dic_engine R.SU / R.SV): the Gauss-Newton
% covariance must predict the ACTUAL scatter of the displacement error.
% Synthetic speckle, known subpixel shift, independent Gaussian noise of
% 1 % of the intensity range in BOTH frames (an ordinary 8-bit camera):
%   (a) predicted rms sigma vs measured std of the error, U and V;
%   (b) coverage: |error| < 2 sigma at ~95 % of the nodes (Gaussian 95.4);
%   (c) noise-free frames give a sigma far below the noisy one - the
%       estimate follows the noise, not the subset size alone.
% At 2-4 % noise the estimate is known to run 15-20 % low (gradient noise
% inflates the Hessian); that regime is documented, not asserted here.
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(420), fspecial_gauss(2.0), 'same');
    img = img(11:410, 11:410);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 0.30; TV = -0.20;
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    P = struct('subset_size', 41, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.2, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3);
    sn = 0.01;
    R = dic_engine(img + sn * randn(size(img)), shifted + sn * randn(size(img)), P);
    m = isfinite(R.U) & R.CONV;
    eu = R.U(m) - TU; ev = R.V(m) - TV;
    ru = sqrt(mean(R.SU(m).^2)) / std(eu);
    rv = sqrt(mean(R.SV(m).^2)) / std(ev);
    z = [eu ./ R.SU(m); ev ./ R.SV(m)];
    cov2 = 100 * mean(abs(z) < 2);
    R0 = dic_engine(img, shifted, P);
    m0 = isfinite(R0.U) & R0.CONV;
    q0 = sqrt(mean(R0.SU(m0).^2)) / sqrt(mean(R.SU(m).^2));

    rows = {};
    rows = chk(rows, 'nodes used', nnz(m), 900, Inf, '%.0f');
    rows = chk(rows, 'predicted / measured sigma_u at 1 % noise', ru, 0.85, 1.15, '%.3f');
    rows = chk(rows, 'predicted / measured sigma_v at 1 % noise', rv, 0.85, 1.15, '%.3f');
    rows = chk(rows, 'coverage |err| < 2 sigma, %', cov2, 90, 99, '%.1f');
    rows = chk(rows, 'noise-free sigma / noisy sigma', q0, 0, 0.1, '%.3f');
end

function rows = test_multiseed()
% Multi-seed RG propagation (dic_engine P.multi_seed). A black band across
% the FULL frame width - wider than the subset window, so a whole grid row
% is excluded by the black-region mask - splits the speckle grid into two
% DISCONNECTED regions: the synthetic version of a weld bead / tape /
% glare barrier. The classical single front can never reach the far side
% (it stays NaN, as on the NLMK weld images,
% 2026-07-14); the re-seed restart must recover it with the same rigid
% shift. The seed is pinned to the TOP region (P.seed_xy) so "near" and
% "far" sides are deterministic. Checks:
%   (a) single-seed run (multi_seed=false) leaves the far side empty,
%       i.e. the fixture is a real barrier;
%   (b) multi-seed run recovers BOTH sides with the correct displacement;
%   (c) on barrier-free data multi-seed is bit-identical to single-seed
%       (no re-seed may trigger when one front reaches everything).
% Deterministic (rng 0), toolbox-free, no external data.
    addpath(proj('6_DIC'));
    rng(0);
    img = conv2(rand(320), fspecial_gauss(2.0), 'same');
    img = img(11:310, 11:310);
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    img = 0.2 + 0.8 * img;              % keep speckle above black_thr
    barr = img;  barr(131:170, :) = 0;  % 40-px black band, full width:
    % subset 21 px -> grid nodes at rows 141..160 have an all-black window
    % (movmax < black_thr) and are masked out, cutting the grid in two.
    TU = 3.30; TV = -2.20;
    Fi = griddedInterpolant({1:size(barr,1), 1:size(barr,2)}, barr, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(barr,2), 1:size(barr,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(barr));
    P = struct('subset_size', 21, 'grid_step', 10, 'search_range', 15, ...
        'black_thr', 0.15, 'zncc_min', 0.5, 'max_iter', 30, 'conv_tol', 1e-3, ...
        'seed_xy', [150 40]);           % pin the seed to the TOP region

    % (a) classical single seed: far side of the band must stay NaN
    P1 = P; P1.multi_seed = false;
    R1 = dic_engine(barr, shifted, P1);
    far1 = R1.U(R1.gy(:) > 170, :);          % nodes clearly beyond the band
    near1 = R1.U(R1.gy(:) < 130, :);
    blocked = ~any(isfinite(far1(:))) && any(isfinite(near1(:)));

    % (b) multi-seed (default): both sides recovered, correct shift
    R2 = dic_engine(barr, shifted, P);
    far2 = R2.U(R2.gy(:) > 170, :);
    fv = far2(isfinite(far2));
    eu = mean(R2.U(:), 'omitnan') - TU;
    ev = mean(R2.V(:), 'omitnan') - TV;

    % (c) barrier-free frame: multi-seed must change nothing (no restart)
    Fi0 = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    sh0 = reshape(Fi0(YY(:) - TV, XX(:) - TU), size(img));
    Ra = dic_engine(img, sh0, P1);           % single-seed
    Rb = dic_engine(img, sh0, P);            % multi-seed allowed
    inert = isequaln(Ra.U, Rb.U) && isequaln(Ra.V, Rb.V) && Rb.nseeds == 1;

    rows = {};
    rows = chk(rows, 'single seed: far side blocked (fixture valid)', double(blocked), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'multi-seed: seeds used', R2.nseeds, 2, 8, '%.0f');
    rows = chk(rows, 'multi-seed: far-side nodes recovered', numel(fv), 100, Inf, '%.0f');
    rows = chk(rows, 'multi-seed: far-side |u err|, px', abs(mean(fv) - TU), 0, 0.05, '%.5f');
    rows = chk(rows, 'multi-seed: full-field |err u|, px', abs(eu), 0, 0.05, '%.5f');
    rows = chk(rows, 'multi-seed: full-field |err v|, px', abs(ev), 0, 0.05, '%.5f');
    rows = chk(rows, 'barrier-free: multi-seed inert (bit-for-bit)', double(inert), 0.5, 1.5, '%.0f');
end

function rows = test_extensometer()
% Virtual extensometer (virtual_extensometer.m): gauge-length readings
% across arbitrary point pairs - the "number in the middle" for a zone DIC
% cannot track directly (weld bead / crack / bond line), read from speckle
% on the two banks alone. A pure post-processor unit test on ANALYTIC
% displacement fields with exactly known answers: no images, deterministic,
% toolbox-free.
%   (a) rigid-body translation -> zero strain / zero opening on any gauge;
%   (b) uniform x-stretch eps=0.01 -> an x-gauge reads eps and the relative
%       du exactly, a y-gauge in the SAME field reads ~0 (orientation is
%       honoured - opening/shear are per-gauge-axis);
%   (c) a pure separation D across a horizontal seam -> a gauge straddling
%       the seam (one end on each bank, the seam itself never sampled) reads
%       opening = D, returned in mm under a scale;
%   (d) a NaN at a gauge end is INVALID with snapping off and RECOVERED by
%       the default nearest-valid-node snap;
%   (e) a 2-state stack returns a 2-column reading curve.
    addpath(proj('6_DIC'));
    gx = 10:10:300;  gy = 10:10:300;
    [X, Y] = meshgrid(gx, gy);
    Z = zeros(size(X));

    % (a) rigid body: U=+5, V=-3 everywhere; horizontal + vertical gauge
    Ea = virtual_extensometer(gx, gy, Z+5, Z-3, ...
        [50 150 250 150; 150 50 150 250]);
    rb_strain = max(abs(Ea.strain));
    rb_open   = max(abs(Ea.opening));

    % (b) uniform x-stretch eps = 0.01: U = eps*X, V = 0
    epsx = 0.01;
    Eb = virtual_extensometer(gx, gy, epsx*X, Z, ...
        [50 150 250 150; ...      % x-gauge, L0 = 200 px
         150 50 150 250]);        % y-gauge in the same field -> ~0
    ex_strain = Eb.strain(1);
    ex_du     = Eb.du(1);         % expect eps*(250-50) = 2.0 px
    ey_strain = abs(Eb.strain(2));

    % (c) pure separation D across a horizontal seam at Y = 150
    D = 4;  mmpx = 0.02;
    Vc = Z;  Vc(Y > 150) = D;     % bottom bank displaced +D in y
    Ec = virtual_extensometer(gx, gy, Z, Vc, [150 100 150 200], 'Scale', mmpx);
    open_mm = Ec.opening(1);      % expect D*mmpx = 0.08 mm
    open_px = Ec.opening(1)/mmpx; % expect D = 4 px

    % (d) NaN at a gauge end: invalid with snap off, recovered by default
    Ud = epsx*X;  Vd = Z;
    [~, jn] = min(abs(gx - 250));  [~, in] = min(abs(gy - 150));
    Ud(in, jn) = NaN;  Vd(in, jn) = NaN;      % kill the exact end node
    En = virtual_extensometer(gx, gy, Ud, Vd, [50 150 250 150], 'SnapRadius', 0);
    Es = virtual_extensometer(gx, gy, Ud, Vd, [50 150 250 150]);   % default snap
    off_valid = double(En.valid(1));          % expect 0
    on_valid  = double(Es.valid(1));          % expect 1

    % (e) 2-state stack -> 2-column curve (rigid state, then stretch state)
    Ek = virtual_extensometer(gx, gy, cat(3, Z+5, epsx*X), cat(3, Z-3, Z), ...
        [50 150 250 150]);
    ncol   = size(Ek.strain, 2);
    stk_ok = abs(Ek.strain(1,1)) < 1e-9 && abs(Ek.strain(1,2) - epsx) < 1e-6;

    rows = {};
    rows = chk(rows, '(a) rigid body: max|strain|',        rb_strain, 0, 1e-9, '%.2e');
    rows = chk(rows, '(a) rigid body: max|opening| px',    rb_open,   0, 1e-9, '%.2e');
    rows = chk(rows, '(b) x-stretch: axial strain',        ex_strain, 0.00999, 0.01001, '%.6f');
    rows = chk(rows, '(b) x-stretch: relative du px',      ex_du,     1.999, 2.001, '%.4f');
    rows = chk(rows, '(b) y-gauge in x-field: |strain|',   ey_strain, 0, 1e-9, '%.2e');
    rows = chk(rows, '(c) seam opening, mm',               open_mm,   0.0799, 0.0801, '%.4f');
    rows = chk(rows, '(c) seam opening, px (= D)',         open_px,   3.999, 4.001, '%.4f');
    rows = chk(rows, '(d) NaN end invalid (snap off)',     off_valid, -0.5, 0.5, '%.0f');
    rows = chk(rows, '(d) NaN end recovered (snap on)',    on_valid,  0.5, 1.5, '%.0f');
    rows = chk(rows, '(e) stack yields 2-col curve',       ncol,      2, 2, '%.0f');
    rows = chk(rows, '(e) stack readings correct',         double(stk_ok), 0.5, 1.5, '%.0f');
end

function rows = test_dice_gold()
% Validation against the bundled official DICe gold result (frame 01).
% Regression intervals around the values measured on R2025b, 2026-06-12:
% bias u +0.0077, bias v -0.0049, RMS u 0.0086, RMS v 0.0064 px.
    evalc('run(proj(''9_Validation_DICe'', ''validate_dice.m''))');
    rows = {};
    rows = chk(rows, 'gold frame compared', double(strcmp(stats(1).label, 'official gold')), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'common points', stats(1).npts, 900, 1000, '%.0f');
    rows = chk(rows, 'bias u vs DICe, px', stats(1).bias_u, 0.004, 0.012, '%.4f');
    rows = chk(rows, 'bias v vs DICe, px', stats(1).bias_v, -0.009, -0.001, '%.4f');
    rows = chk(rows, 'RMS u vs DICe, px', stats(1).rms_u, 0.005, 0.013, '%.4f');
    rows = chk(rows, 'RMS v vs DICe, px', stats(1).rms_v, 0.003, 0.010, '%.4f');
    rows = chk(rows, 'all compared frames pass', double(all([stats.pass])), 0.5, 1.5, '%.0f');
end

function rows = test_stereo_cal()
% Toolbox-free camera calibration (11_Stereo, phase-2 stereo groundwork).
% A fast cam0 subset (8 spread views) recomputes intrinsics from the
% tracked AMCalB dot-target images: detect_dots -> index_grid ->
% calibrate_zhang. The pipeline is deterministic (no rng), so the
% intervals are pinned around the values first measured on R2025b,
% 2026-07-05: fx 6693.2, fy/fx 1.0001, u0 957.6, v0 655.8, RMS 0.315 px,
% 136 dots + 3 donuts on view 0. The recovered focal must also agree with
% DICe's own rig calibration (cal.xml fx 6636) to within a few percent,
% which is what confirms AMCalB is the GT4 stereo rig.
    addpath(proj('11_Stereo'));
    caldir = proj('9_Validation_DICe', 'calibration');
    D = detect_dots(imread(fullfile(caldir, 'AMCalB-0000_0.jpeg')));
    vs = [0 8 16 24 32 40 48 56];
    R = calibrate_cameras(struct('cams', 0, 'views', vs, 'save', false, ...
        'verbose', false, 'caldir', caldir));
    C = R.cam(1);
    fx = C.K(1,1); fy = C.K(2,2); u0 = C.K(1,3); v0 = C.K(2,3);
    rows = {};
    rows = chk(rows, 'view0: donut fiducials', size(D.donut_xy,1), 3, 3, '%.0f');
    rows = chk(rows, 'view0: dots detected', D.n, 120, 160, '%.0f');
    rows = chk(rows, 'subset views calibrated', numel(C.kept), 8, 8, '%.0f');
    rows = chk(rows, 'fx, px', fx, 6650, 6735, '%.1f');
    rows = chk(rows, 'fy/fx ratio', fy/fx, 0.98, 1.02, '%.4f');
    rows = chk(rows, 'u0 in central 20%, px', u0, 768, 1152, '%.1f');
    rows = chk(rows, 'v0 in central 20%, px', v0, 480, 720, '%.1f');
    rows = chk(rows, 'reprojection RMS, px', C.rms, 0, 0.45, '%.4f');
    rows = chk(rows, 'fx agrees with DICe cal.xml (<3%)', abs(fx-6636)/6636, 0, 0.03, '%.4f');
end

function rows = test_stereo3()
% Stereo extrinsics + triangulation (11_Stereo, phase-2 groundwork). A small
% both-camera calibration (9 spread views) yields the relative pose
% cam0->cam1; the view-0 grid is triangulated and checked for planarity and
% uniform reconstructed pitch. Anchoring the (unit-pitch) baseline to DICe's
% cal.xml mm baseline recovers the physical dot pitch. Deterministic (no
% rng); intervals pinned to the R2025b, 2026-07-05 measured run: rel
% rotation 18.71 deg, baseline CV 0.54%, planarity 0.0067 (0.7% of a pitch),
% reconstructed pitch 0.999, dot pitch 7.05 mm (which refutes the 1.90425 mm
% EXIF guess; the full 66-view solution gives 6.97 mm).
    addpath(proj('11_Stereo'));
    caldir = proj('9_Validation_DICe', 'calibration');
    vs = [0 8 16 24 32 40 48 56 64];
    R = calibrate_cameras(struct('cams', [0 1], 'views', vs, 'save', false, ...
        'verbose', false, 'caldir', caldir));
    E = stereo_extrinsics(R);
    D0 = detect_dots(imread(fullfile(caldir, 'AMCalB-0000_0.jpeg')));
    D1 = detect_dots(imread(fullfile(caldir, 'AMCalB-0000_1.jpeg')));
    G0 = index_grid(D0); G1 = index_grid(D1);
    k0 = G0.i*1000 + G0.j; k1 = G1.i*1000 + G1.j;
    [~, ia, ib] = intersect(k0, k1);
    X = triangulate_dlt(G0.uv(ia,:), G1.uv(ib,:), R.K0, R.dist0, R.K1, R.dist1, E.R, E.T);
    c = mean(X,1); Xc = X - c; [~,~,Vp] = svd(Xc,0);
    planrms = sqrt(mean((Xc*Vp(:,3)).^2));
    ii = G0.i(ia); jj = G0.j(ia); key = k0(ia); dd = [];
    for q = 1:numel(ii)
        for de = [1 0; 0 1].'
            t = find(key == (ii(q)+de(1))*1000 + (jj(q)+de(2)), 1);
            if ~isempty(t), dd(end+1) = norm(X(q,:) - X(t,:)); end %#ok<AGROW>
        end
    end
    pitch_recon = mean(dd); pitch_cv = std(dd)/mean(dd);
    base_mm = dice_stereo_baseline(proj('9_Validation_DICe', 'stereo_d_spec', 'cal.xml'));
    pitch_mm = base_mm / E.baseline;
    rows = {};
    rows = chk(rows, 'common stereo views', E.nviews, 8, 9, '%.0f');
    rows = chk(rows, 'rel rotation cam0->cam1, deg', E.ang_deg, 18.3, 19.1, '%.3f');
    rows = chk(rows, 'baseline consistency CV, %', E.T_cv*100, 0, 1.5, '%.3f');
    rows = chk(rows, 'triangulation planarity RMS, pitch', planrms, 0, 0.02, '%.4f');
    rows = chk(rows, 'reconstructed pitch, units', pitch_recon, 0.99, 1.01, '%.4f');
    rows = chk(rows, 'reconstructed pitch CV, %', pitch_cv*100, 0, 1.0, '%.3f');
    rows = chk(rows, 'anchored dot pitch (DICe mm), mm', pitch_mm, 6.5, 7.6, '%.3f');
end

function rows = test_stereo_gold()
% Stereo (3D) DIC vs DICe's stereo gold (11_Stereo, phase-2). The same
% deterministic 9-view both-camera calibration as test_stereo3 builds the
% rig geometry in-memory (nothing saved); STEREO_DIC then processes GT4
% frames 0 and 5 of stereo_d_spec and is compared against DICe's own gold
% solutions at the gold subset centres. Intervals pinned to this test's
% first measured run (R2025b, 2026-07-06); the full-view geometry numbers
% (validate_dice_stereo, same date) are better still - e.g. gold-pairs
% epipolar 0.071 px vs 0.775 with the 9-view geometry here. NOTE: gold
% CROSS_CORR is frame-constant (DICe cross-matches the cameras once), so
% disparity is only apples-to-apples at frame 0; frame 5 checks the
% temporal channel.
    addpath(proj('11_Stereo'), proj('6_DIC'));
    caldir = proj('9_Validation_DICe', 'calibration');
    spec   = proj('9_Validation_DICe', 'stereo_d_spec');
    vs = [0 8 16 24 32 40 48 56 64];
    R = calibrate_cameras(struct('cams', [0 1], 'views', vs, 'save', false, ...
        'verbose', false, 'caldir', caldir));
    G = stereo_geometry(struct('cal', R, 'dicecal', fullfile(spec, 'cal.xml')));

    P = struct('subset_size', 25, 'grid_step', 40, 'search_range', 40, ...
        'search_temporal', 10, 'black_thr', 0.15, 'zncc_min', 0.75, ...
        'zncc_min_cross', 0.50, 'max_iter', 30, 'conv_tol', 1e-3);
    f0 = {fullfile(spec, 'images', 'GT4-0000_0.tif'), ...
          fullfile(spec, 'images', 'GT4-0005_0.tif')};
    f1 = {fullfile(spec, 'images', 'GT4-0000_1.tif'), ...
          fullfile(spec, 'images', 'GT4-0005_1.tif')};
    S = stereo_dic(f0, f1, G, P);

    % ---- frame 0: disparity + geometry + 3D ----
    T0 = read_dice_gold(fullfile(spec, 'gold', 'DICe_solution_0.txt'));
    v = T0.MATCH >= 0;
    xg = T0.COORDINATE_X(v); yg = T0.COORDINATE_Y(v);
    qg = T0.CROSS_CORR_Q(v); rg = T0.CROSS_CORR_R(v);
    D0 = S.frames(1).D;
    Fu = griddedInterpolant({D0.gy, D0.gx}, D0.U, 'linear', 'none');
    Fv = griddedInterpolant({D0.gy, D0.gx}, D0.V, 'linear', 'none');
    ql = Fu(yg, xg); rl = Fv(yg, xg);
    m = isfinite(ql) & isfinite(rl);
    dq = ql(m) - qg(m); dr = rl(m) - rg(m);
    med0 = median(hypot(dq, dr));
    % Epipolar residual of the gold pairs under TWO geometries. Under
    % DICe's own cal.xml (parsed, relative pose via dice_stereo_baseline)
    % it checks the parsing + Q/R convention with a precise reference
    % geometry. Under our 9-view subset geometry it measures THAT
    % geometry's accuracy floor (first measured 0.775 px - a ~0.007 deg
    % rotation error at f~6640 px; the full 66-view calibration measured
    % 0.071 px on the same quantity, R2025b, 2026-07-06).
    epi_gold = epi_rms_local([xg(m) yg(m)], [xg(m)+qg(m) yg(m)+rg(m)], G.F);
    Cd = read_dice_cal(fullfile(spec, 'cal.xml'));
    [~, Rd, Td] = dice_stereo_baseline(Cd);
    Txd = [0 -Td(3) Td(2); Td(3) 0 -Td(1); -Td(2) Td(1) 0];
    Fd = (Cd.cam(2).K.') \ (Txd * Rd) / Cd.cam(1).K;
    epi_dice = epi_rms_local([xg(m) yg(m)], [xg(m)+qg(m) yg(m)+rg(m)], Fd);
    Xl = triangulate_dlt([xg(m) yg(m)], [xg(m)+ql(m) yg(m)+rl(m)], ...
        G.K0, G.d0, G.K1, G.d1, G.R, G.T) * G.pitch_mm;
    Xd = triangulate_dlt([xg(m) yg(m)], [xg(m)+qg(m) yg(m)+rg(m)], ...
        G.K0, G.d0, G.K1, G.d1, G.R, G.T) * G.pitch_mm;
    rms3d = sqrt(mean(sum((Xl - Xd).^2, 2)));

    % ---- frame 5: temporal channel ----
    T5 = read_dice_gold(fullfile(spec, 'gold', 'DICe_solution_5.txt'));
    v5 = T5.MATCH >= 0;
    x5 = T5.COORDINATE_X(v5); y5 = T5.COORDINATE_Y(v5);
    Tu = griddedInterpolant({S.gy, S.gx}, S.frames(2).U0, 'linear', 'none');
    Tv = griddedInterpolant({S.gy, S.gx}, S.frames(2).V0, 'linear', 'none');
    ul = Tu(y5, x5); vl = Tv(y5, x5);
    m5 = isfinite(ul) & isfinite(vl);
    ug = T5.DISPLACEMENT_X(v5); vg = T5.DISPLACEMENT_Y(v5);
    du = ul(m5) - ug(m5); dv = vl(m5) - vg(m5);
    rms_u = sqrt(mean(du.^2)); rms_v = sqrt(mean(dv.^2));

    % Intervals pinned to the first measured run of this test (R2025b,
    % 2026-07-06: pts 103, epi9 0.775, median 0.050, RMS 0.103/0.059,
    % temporal 0.0073/0.0059, 3D 0.0190 mm) with ~50% headroom; the DICe-
    % geometry epipolar row is provisional (0..0.5) until first measured.
    rows = {};
    rows = chk(rows, 'matched points at gold centres', nnz(m), 95, 115, '%.0f');
    rows = chk(rows, 'gold pairs epipolar RMS (DICe cal.xml geometry), px', epi_dice, 0, 0.5, '%.3f');
    rows = chk(rows, 'gold pairs epipolar RMS (9-view geometry), px', epi_gold, 0, 1.5, '%.3f');
    rows = chk(rows, 'frame0 disparity vs gold: median, px', med0, 0, 0.10, '%.3f');
    rows = chk(rows, 'frame0 disparity vs gold: RMS x, px', sqrt(mean(dq.^2)), 0, 0.15, '%.3f');
    rows = chk(rows, 'frame0 disparity vs gold: RMS y, px', sqrt(mean(dr.^2)), 0, 0.12, '%.3f');
    rows = chk(rows, 'frame5 temporal vs gold: RMS x, px', rms_u, 0, 0.015, '%.4f');
    rows = chk(rows, 'frame5 temporal vs gold: RMS y, px', rms_v, 0, 0.015, '%.4f');
    rows = chk(rows, 'frame0 3D surface RMS vs gold, mm', rms3d, 0, 0.04, '%.4f');
end

function r = epi_rms_local(uv0, uv1, F)
% RMS distance of uv1 to the epipolar lines F*uv0 (px).
    n = size(uv0, 1);
    l = F * [uv0, ones(n,1)].';
    num = sum([uv1, ones(n,1)].' .* l, 1);
    den = hypot(l(1,:), l(2,:));
    r = sqrt(mean((num ./ den).^2));
end

function rows = test_surf_strain()
% Surface (membrane) strain from stereo 3D points (surf_strain.m). Purely
% synthetic and deterministic - no images needed. Checks the defining
% properties of the algorithm: (1) an exactly known in-plane affine
% deformation of a tilted plane is recovered to machine precision
% (principal strains, frame-independent); (2) rigid-body motion of a
% curved (cylindrical) surface gives zero strain; (3) a uniform axial
% stretch lambda of the cylinder gives principal strains
% [(lambda^2-1)/2, 0]; (4) all of it survives a NaN hole in the grid.
    rows = {};
    addpath(proj('11_Stereo'));

    % --- 1) tilted flat plane, exact in-plane affine deformation ---
    [gx, gy] = meshgrid((0:39)*1.3, (0:29)*1.7);
    R0 = rotv([1; 2; 0.5], 0.4);
    A  = [gx(:), gy(:), zeros(numel(gx),1)] * R0.';   % plane in 3D
    e1g = R0(:,1); e2g = R0(:,2); ng = R0(:,3);
    Fin = [1.013 0.004; 0.004 0.9945];                % in-plane gradient
    loc  = [A*e1g, A*e2g];
    locd = loc * Fin.';
    Ad = locd(:,1)*e1g.' + locd(:,2)*e2g.' + (A*ng)*ng.';
    Xref = reshape(A,  [size(gx) 3]);
    Xdef = reshape(Ad, [size(gx) 3]);
    E = surf_strain(Xref, Xdef);
    lam_gt = sort(eig(0.5*(Fin.'*Fin - eye(2))));
    err = 0;
    for i = 4:size(gx,1)-3
        for j = 4:size(gx,2)-3
            lam = sort(eig([E.exx(i,j) E.exy(i,j); E.exy(i,j) E.eyy(i,j)]));
            err = max(err, max(abs(lam - lam_gt)));
        end
    end
    rows = chk(rows, 'flat affine: principal-strain error', err, 0, 1e-9, '%.1e');

    % --- 2) cylinder patch, rigid-body motion -> zero strain ---
    [TH, Z] = meshgrid(linspace(-0.5, 0.5, 45), linspace(0, 60, 35));
    Rc = 80;
    pts = cat(3, Rc*sin(TH), Z, Rc*(1-cos(TH)));
    Rb = rotv([0.2; 1; 0.1], 0.25);
    Xd = reshape(reshape(pts, [], 3) * Rb.' + [5 -3 2], size(pts));
    E = surf_strain(pts, Xd);
    mx = max(abs([E.exx(:); E.eyy(:); E.exy(:)]), [], 'omitnan');
    rows = chk(rows, 'cylinder rigid body: max |strain|', mx, 0, 1e-8, '%.1e');

    % --- 3) cylinder axial stretch: principal strains vs analytic ---
    lam = 1.008;
    Xd = pts; Xd(:,:,2) = Xd(:,:,2) * lam;
    E = surf_strain(pts, Xd);
    E1gt = 0.5*(lam^2 - 1);
    err = 0;
    for i = 4:size(pts,1)-3
        for j = 4:size(pts,2)-3
            v = sort(eig([E.exx(i,j) E.exy(i,j); E.exy(i,j) E.eyy(i,j)]));
            err = max(err, max(abs(v - [0; E1gt])));
        end
    end
    rows = chk(rows, 'cylinder stretch: error vs (lam^2-1)/2', err, 0, 1e-6, '%.1e');

    % --- 4) NaN hole in the grid: still correct around it ---
    ptsh = pts; Xdh = Xd;
    ptsh(11:14, 21:26, :) = NaN; Xdh(11:14, 21:26, :) = NaN;
    E = surf_strain(ptsh, Xdh);
    err = 0; nv = 0;
    for i = 4:size(pts,1)-3
        for j = 4:size(pts,2)-3
            if ~isfinite(E.exx(i,j)), continue; end
            nv = nv + 1;
            v = sort(eig([E.exx(i,j) E.exy(i,j); E.exy(i,j) E.eyy(i,j)]));
            err = max(err, max(abs(v - [0; E1gt])));
        end
    end
    rows = chk(rows, 'with a hole: strain error', err, 0, 1e-6, '%.1e');
    rows = chk(rows, 'with a hole: valid interior nodes', nv, 800, 1e9, '%.0f');
end

function R = rotv(ax, ang)
% Rotation matrix about axis AX by ANG rad (Rodrigues).
    ax = ax(:) / norm(ax);
    K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
    R = eye(3) + sin(ang)*K + (1 - cos(ang))*(K*K);
end

function rows = test_gt_gen()
% Synthetic ground-truth factory (12_Datasets/gen_dataset.m). The imaging
% chain is interpolation-free by construction (analytic Gaussian speckles,
% exact inverse warp), so the checks are sharp: the reported inverse-warp
% residual must be at machine level, the untouched 2D engine must recover
% a noise-free affine field to a few hundredths of a pixel (first-order
% subset shape function matches an affine field exactly), and the same
% seed must reproduce the images bit-for-bit.
    rows = {};
    addpath(proj('12_Datasets'), proj('6_DIC'));
    od = tempname;
    M = gen_dataset(struct('outdir', od, 'n_pairs', 2, 'imsize', [256 256], ...
        'field', 'affine', 'noise_std', 0, 'illum_grad', 0, 'seed', 42));
    rows = chk(rows, 'inverse-warp residual, px', ...
        max([M.qa.inv_residual]), 0, 1e-8, '%.1e');

    I0 = double(imread(fullfile(od, 'pair_0001_ref.png'))) / 65535;
    I1 = double(imread(fullfile(od, 'pair_0001_def.png'))) / 65535;
    S  = load(fullfile(od, 'pair_0001_gt.mat'));
    P = struct('subset_size', 21, 'grid_step', 12, 'search_range', 8, ...
        'black_thr', 0.05, 'zncc_min', 0.6, 'max_iter', 30, 'conv_tol', 1e-3);
    R = dic_engine(I0, I1, P);
    [GX, GY] = meshgrid(R.gx, R.gy);
    idx = sub2ind(size(S.U), GY, GX);      % GT is 0-based-coord array: U(gy,gx)
    m = ~isnan(R.U);
    du = R.U(m) - S.U(idx(m));
    dv = R.V(m) - S.V(idx(m));
    rows = chk(rows, 'engine vs GT: points', nnz(m), 100, 1e9, '%.0f');
    rows = chk(rows, 'engine vs GT: RMS u, px', sqrt(mean(du.^2)), 0, 0.02, '%.4f');
    rows = chk(rows, 'engine vs GT: RMS v, px', sqrt(mean(dv.^2)), 0, 0.02, '%.4f');
    rows = chk(rows, 'engine vs GT: |bias u|, px', abs(mean(du)), 0, 0.01, '%.4f');

    od2 = tempname;
    gen_dataset(struct('outdir', od2, 'n_pairs', 1, 'imsize', [256 256], ...
        'field', 'affine', 'noise_std', 0, 'illum_grad', 0, 'seed', 42));
    same = isequal(imread(fullfile(od,  'pair_0001_ref.png')), ...
                   imread(fullfile(od2, 'pair_0001_ref.png'))) && ...
           isequal(imread(fullfile(od,  'pair_0001_def.png')), ...
                   imread(fullfile(od2, 'pair_0001_def.png')));
    rows = chk(rows, 'same seed => bit-identical images', double(same), 0.5, 1.5, '%.0f');
    try rmdir(od, 's');  catch, end
    try rmdir(od2, 's'); catch, end
end

function rows = test_validator()
% External-field validator (13_Validator). Round trip: the DICe gold file
% itself, read by the universal importer and validated against itself,
% must show zero error and PASS - this exercises the comment skipping,
% header-name mapping, SIGMA filtering, interpolation and verdict logic.
% A header-less shuffled CSV of the same points must give the same PASS.
    rows = {};
    addpath(proj('13_Validator'));
    gold = proj('9_Validation_DICe', 'dice_reference', 'DICe_solution_0.txt');
    E = import_external_dic(gold);
    rows = chk(rows, 'gold parsed: points', E.info.n_points, 900, 1100, '%.0f');

    od = tempname; mkdir(od);
    R = validate_external(gold, struct('outdir', od));
    rows = chk(rows, 'round-trip |bias u|, px', abs(R.bias_u), 0, 1e-9, '%.1e');
    rows = chk(rows, 'round-trip RMS v, px', R.rms_v, 0, 1e-9, '%.1e');
    rows = chk(rows, 'round-trip verdict = PASS', double(R.pass), 0.5, 1.5, '%.0f');

    rng(1, 'twister');
    T = [E.x, E.y, E.u, E.v];
    T = T(randperm(size(T, 1)), :);
    f2 = fullfile(od, 'plain.csv');
    fid = fopen(f2, 'w');
    fprintf(fid, '%.9g,%.9g,%.9g,%.9g\n', T.');
    fclose(fid);
    R2 = validate_external(f2, struct('outdir', od));
    rows = chk(rows, 'header-less CSV verdict = PASS', double(R2.pass), 0.5, 1.5, '%.0f');
    try rmdir(od, 's'); catch, end
end

function rows = test_stereo_dataset(dataDir)
% OPT-IN welded-plate stereo regression (RUN_TESTS('stereo-data')). Runs the
% full validation END-TO-END on the local dataset and freezes the numbers
% from the 2026-07-14 validation (multi-seed engine, speckle-zone verdict).
% NOT in the base suite: it depends on an external dataset and takes minutes.
% Layouts supported: images mixed in dataDir (historical, auto-split), or
% the 2026-07-13 relayout with Dataset\ + Calibration\ subfolders (21r.csv
% at the dataset root is copied next to the results if needed).
% Frozen reference - our handheld Zhang calibration at the true 7 mm pitch,
% vs the Vic-3D L0 line export (21r.csv). Re-frozen 2026-08-13 after the
% DETECT_STRIP fix, ROI now covers both speckled zones (rows 1..3856,
% was 1..1804):
%   cam0 reprojection RMS 0.1546 px, cam1 0.1631 px, baseline 274.31 mm;
%   SPECKLE-ZONE depth agreement s_z 0.9917 with worst per-state rms 0.6 um
%   over 123 rows spanning 99% of the reference line
%   (the all-row s_z is NOT gated: over the bare weld zone Vic's export is
%   an interpolated bridge, not a measurement - see compare_vic3d header);
%   multi-seed keeps >= 4300 valid nodes per deformed state (was ~2800
%   when the single front died at the weld line).
% The calibration residuals also moved (0.2170/0.2234 -> 0.1546/0.1631), and
% not because of the ROI (not used in calibration): they were frozen before
% the 2026-07-30 dot-centre change and not re-measured after it.
% Coverage is checked separately, because every other number in
% this test is computed on the overlap between the two lines and therefore
% cannot see a specimen that was only half measured.
% Intervals carry modest headroom around the measured values.
    addpath(proj('11_Stereo'), proj('6_DIC'));
    if ~isfolder(dataDir)
        error(['stereo dataset not found at %s - move it there or call ' ...
               'RUN_TESTS(''stereo-data'', ''<path to DIC 3D dataset>'').'], dataDir);
    end
    if ~isfile(fullfile(dataDir, '21r.csv'))
        error(['%s has no 21r.csv (the Vic-3D L0 export) - nothing to ' ...
               'compare against.'], dataDir);
    end

    % layout: split Dataset\/Calibration\ subfolders, or mixed single folder
    if isfolder(fullfile(dataDir, 'Dataset')) && ...
            isfolder(fullfile(dataDir, 'Calibration'))
        runDir = fullfile(dataDir, 'Dataset');
        calDir = fullfile(dataDir, 'Calibration');
        if ~isfile(fullfile(runDir, '21r.csv'))
            copyfile(fullfile(dataDir, '21r.csv'), runDir);
        end
        if isfile(fullfile(dataDir, 'project.xml')) && ...
                ~isfile(fullfile(runDir, 'project.xml'))
            copyfile(fullfile(dataDir, 'project.xml'), runDir);
        end
        evalc('run_stereo_dataset(runDir, ''Pitch'', 7, ''CalDir'', calDir)');
    else
        runDir = dataDir;
        evalc('run_stereo_dataset(dataDir, ''Pitch'', 7)');
    end
    outDir = fullfile(runDir, 'litedic_out');
    out = compare_vic3d(runDir);

    cal = load(fullfile(outDir, 'cal_stereo.mat'), 'rms0', 'rms1');
    Lr  = load(fullfile(outDir, 'stereo_results.mat'), 'G', 'S');
    baseline_mm = norm(Lr.G.T);      % pitch=7 -> object points in mm -> mm
    nv = arrayfun(@(f) f.nvalid, Lr.S.frames(2:end));

    rows = {};
    rows = chk(rows, 'cam0 reprojection RMS, px', cal.rms0, 0.140, 0.175, '%.4f');
    rows = chk(rows, 'cam1 reprojection RMS, px', cal.rms1, 0.148, 0.182, '%.4f');
    rows = chk(rows, 'stereo baseline, mm',       baseline_mm, 271.0, 276.0, '%.2f');
    rows = chk(rows, 'speckle-zone depth agreement s_z', out.s_z_speckle, 0.955, 1.010, '%.3f');
    rows = chk(rows, 'speckle-zone rms (worst state), mm', out.rms_max_speckle, 0, 0.004, '%.4f');
    rows = chk(rows, 'speckle rows used', out.n_speckle, 80, Inf, '%.0f');
    % catches a specimen measured only in part (rows above only see the
    % measured part)
    rows = chk(rows, 'reference line covered', out.cover, 0.90, 1.001, '%.3f');
    rows = chk(rows, 'multi-seed: min valid nodes per state', min(nv), 4000, Inf, '%.0f');
end
