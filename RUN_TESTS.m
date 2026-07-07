function RUN_TESTS(mode)
% RUN_TESTS - automated regression harness: runs every analysis and checks
% the key numbers against reference intervals.
%
% Usage from the project folder:
%       RUN_TESTS            % main set, ~1-2 min
%       RUN_TESTS('full')    % + FE beam (requires the PDE Toolbox)
%
% The report is printed, saved to TEST_REPORT.txt and copied to the
% clipboard.

if nargin < 1, mode = 'base'; end
do_full = strcmpi(mode, 'full');
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

G = run_case(G, 'CLT (plate_clt.m)',            @test_clt);
G = run_case(G, 'FSDT (plate_fsdt.m)',          @test_fsdt);
G = run_case(G, 'Vibration (plate_vibration.m)', @test_wave);
G = run_case(G, 'DDM cracking (matrix_cracking_ddm.m)', @test_crack);
G = run_case(G, 'Shell (cyl_shell.m)',          @test_cyl);
G = run_case(G, 'Rocket (rocket_loads.m)',      @test_rock);
G = run_case(G, 'Sphere (sphere_demo.m)',       @test_sphere);
G = run_case(G, 'DIC coarse (dic_coarse.m)',    @test_corr);
G = run_case(G, 'DIC subpixel (dic_icgn.m)',    @test_dicpro);
G = run_case(G, 'DIC batch kernel (dic_engine_batch.m)', @test_batch);
G = run_case(G, 'DIC bias S-curve (bias_scurve.m)', @test_scurve);
G = run_case(G, 'DIC mm/px scale (dic_scale.m)', @test_scale);
G = run_case(G, 'DIC vs DICe gold (validate_dice.m)', @test_dice_gold);
G = run_case(G, 'Stereo intrinsics (11_Stereo)', @test_stereo_cal);
G = run_case(G, 'Stereo extrinsics+triangulation (11_Stereo)', @test_stereo3);
G = run_case(G, 'Stereo DIC vs DICe gold (11_Stereo)', @test_stereo_gold);
if do_full
    G = run_case(G, 'FE beam (beam_fem.m, PDE Toolbox)', @test_beam);
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
try
    clipboard('copy', rep);
    fprintf('\nReport copied to the CLIPBOARD and saved to TEST_REPORT.txt.\n');
catch
    fprintf('\nReport saved to TEST_REPORT.txt.\n');
end
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
        fullfile('1_Plates_CLT_FSDT','plate_clt.m'), ...
        fullfile('1_Plates_CLT_FSDT','plate_fsdt.m'), ...
        fullfile('2_Damage_crack','matrix_cracking_ddm.m'), ...
        fullfile('3_Cyl_shell','cyl_shell.m'), ...
        fullfile('4_Plate_vibration','plate_vibration.m'), ...
        fullfile('5_Rocket_strength','rocket_loads.m'), ...
        fullfile('6_DIC','dic_coarse.m'), ...
        fullfile('6_DIC','dic_icgn.m'), ...
        fullfile('6_DIC','dic_engine.m'), ...
        fullfile('6_DIC','dic_engine_batch.m'), ...
        fullfile('6_DIC','pls_strain.m'), ...
        fullfile('6_DIC','dic_scale.m'), ...
        fullfile('6_DIC','scale_from_length.m'), ...
        fullfile('6_DIC','plate_ref.tiff'), ...
        fullfile('7_Sphere_geometry','sphere_demo.m'), ...
        fullfile('9_Validation_DICe','validate_dice.m'), ...
        fullfile('9_Validation_DICe','dice_reference','DICe_solution_0.txt'), ...
        fullfile('10_Challenge2','challenge2_benchmark.m'), ...
        fullfile('8_FEM_beam','beam_fem.m'), ...
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

function rows = test_clt()
    evalc('run(proj(''1_Plates_CLT_FSDT'', ''plate_clt.m''))');
    rows = {};
    rows = chk(rows, 'center deflection, mm', w*1000, 18.85, 19.05, '%.2f');
    rows = chk(rows, 'moment M_x, N*m/m', M(1), -470, -455, '%.1f');
    rows = chk(rows, 'symmetric layup: max|B|, N', max(abs(B(:))), 0, 1e-6, '%.1e');
end

function rows = test_fsdt()
    evalc('run(proj(''1_Plates_CLT_FSDT'', ''plate_fsdt.m''))');
    rows = {};
    rows = chk(rows, 'center deflection, mm', w*1000, 18.85, 19.10, '%.2f');
    rows = chk(rows, 'moment M_x, N*m/m', M(1), -470, -455, '%.1f');
    rows = chk(rows, 'shear stiffness A_s(1,1), N/m', A_s(1,1), 2.4e7, 2.6e7, '%.3g');
    rows = chk(rows, 'shear stiffness A_s(2,2), N/m', A_s(2,2), 2.4e7, 2.6e7, '%.3g');
end

function rows = test_wave()
    evalc('run(proj(''4_Plate_vibration'', ''plate_vibration.m''))');
    rows = {};
    rows = chk(rows, 'frequency f1, Hz', f(1), 5.60, 5.80, '%.3f');
    rows = chk(rows, 'beam check value, Hz', f1_beam, 5.55, 5.67, '%.3f');
    rows = chk(rows, 'f1 above beam value (plate is stiffer)', double(f(1) > f1_beam), 0.5, 1.5, '%.0f');
    rows = chk(rows, 'frequency f5, Hz', f(5), 200, 212, '%.1f');
end

function rows = test_crack()
    evalc('run(proj(''2_Damage_crack'', ''matrix_cracking_ddm.m''))');
    rows = {};
    R = 1 - D2_history;
    rows = chk(rows, 'first cracks, % strain', ...
        strain_steps(find(lambda_history > 1e-3, 1))*100, 0.20, 0.32, '%.2f');
    dR = gradient(R, strain_steps); [~, i_inf] = min(dR);
    rows = chk(rows, 'sigmoid inflection, % strain', strain_steps(i_inf)*100, 0.40, 0.52, '%.2f');
    rows = chk(rows, 'shoulder: E2eff/E2 at 0.30 %', interp1(strain_steps, R, 0.0030), 0.95, 1.001, '%.3f');
    rows = chk(rows, 'lambda at the end of the range', lambda_history(end), 0.88, 1.00, '%.3f');
    rows = chk(rows, 'residual stiffness E2eff/E2', R(end), 0.33, 0.37, '%.3f');
    rows = chk(rows, 'lambda monotonicity', double(all(diff(lambda_history) >= -1e-12)), 0.5, 1.5, '%.0f');
end

function rows = test_cyl()
    evalc('run(proj(''3_Cyl_shell'', ''cyl_shell.m''))');
    rows = {};
    rows = chk(rows, 'membrane sigma_theta, MPa', sigma_theta_mem/1e6, 99.5, 100.5, '%.1f');
    rows = chk(rows, 'membrane sigma_x, MPa', sigma_x_mem/1e6, 49.5, 50.5, '%.1f');
    rows = chk(rows, 'membrane deflection w, mm', w_mem*1e3, 0.41, 0.44, '%.3f');
end

function rows = test_rock()
    evalc('run(proj(''5_Rocket_strength'', ''rocket_loads.m''))');
    rows = {};
    rows = chk(rows, 'tank cylinder: sigma_theta, Pa', sigma_theta_fuel, 0.99e7, 1.01e7, '%.3g');
    rows = chk(rows, 'tank cylinder: sigma_x, Pa', sigma_x_fuel, 4.9e6, 5.1e6, '%.3g');
    rows = chk(rows, 'tank bottom: sigma, Pa', sigma_bottom_fuel, 2.4e6, 2.6e6, '%.3g');
    rows = chk(rows, 'bottom type = flat', double(strcmp(bottom_type_fuel, 'flat')), 0.5, 1.5, '%.0f');
end

function rows = test_sphere()
    evalc('run(proj(''7_Sphere_geometry'', ''sphere_demo.m''))');
    rows = {'ran without errors', 1, ''};
end

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
    % Loading series ("series" shadows a toolbox function, copy first).
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

function rows = test_beam()
    evalc('run(proj(''8_FEM_beam'', ''beam_fem.m''))');
    rows = {'ran without errors (check ~420 Hz on the tip plot)', 1, ''};
end
