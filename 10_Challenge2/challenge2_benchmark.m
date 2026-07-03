% DIC Challenge 2.0 benchmark: measurement resolution, spatial resolution
% and the Metrological Efficiency Indicator (MEI) for LiteDIC.
%
% Method (Reu et al. 2022, "DIC Challenge 2.0", Exp. Mech.):
%   * the Star images encode a VERTICAL sinusoidal field whose period
%     grows linearly across the image:  lambda(x) = 10 + (300-10)/4000
%     * (x-1) px  (Eq. 1). Checked on the images by optical
%     flow: Star 5 is v(x,y) = 0.5*sin(2*pi*y/lambda(x) +
%     phi(x)) px, Star 6 is the strain analogue with v amplitude
%     0.05*lambda(x)/(2*pi), i.e. eps_yy amplitude 0.05. The phase
%     reference phi(x) is x-dependent (a tilted fan), so the analysis
%     below is built to not need it;
%   * measurement resolution n = 1 sigma of the QOI along the center-row
%     line-cut of the UNDEFORMED noise-floor pair (pure noise);
%   * spatial resolution ell_10% = the period at which the measured
%     amplitude first drops 10 % below ground truth, read from a
%     12th-order polynomial fit of amplitude ratio vs period (the raw
%     curve is jagged from Pattern-Induced Bias; the fit is part of the
%     published method);
%   * MEI = n^o * ell_10%, o = 1 for displacement, o = 2 for strain
%     (lower is better). NOTE: the constancy of MEI across
%     the parameter sweep holds for DISPLACEMENT only. The strain MEI
%     falls monotonically with total smoothing (~25x across the published
%     Local.Affine.B subset sweep, ~100x across the subset x VSG sweep),
%     so "mean of 3 lowest" rewards sweep extent; cross-code strain
%     comparisons should be made at matched ell_10% (see the paper,
%     sec-challenge2).
%
% Amplitude extraction: at every 1-px column x the measured field over a
% band of rows is fitted with [sin(2*pi*y/lambda(x)), cos(.), 1] by least
% squares; the first-harmonic amplitude divided by the ground-truth
% amplitude is the attenuation at period lambda(x). This measures exactly
% what the Challenge's per-pixel amplitude data measures (the fit +
% crossing stage is validated against the published Local.Affine.B
% results), is insensitive to the fan's phase reference, and averages
% noise over the band. The noise line-cut itself is evaluated at the
% center row at 1 px spacing straight from the DIC grid - no resampling,
% as required by the Challenge.
%
% The engine only processes a horizontal band of rows around the center
% line; the band contains every row the strain window and the harmonic
% fit need. Points are solved with the batched IC-GN kernel
% (dic_engine_batch.m - same mathematics as dic_engine, all CPU cores,
% no toolboxes); a preflight asserts its agreement with the sequential
% RG-DIC engine on the actual images before anything is measured.
%
% Images: https://sem.org/dicchallenge/ (place the Star folders inside
% "2D-Challenge 2.0" next to this folder - see README.md). If the
% published per-code results CSVs are present, LiteDIC is plotted on top
% of the Local.Affine.B reference curve for direct comparison.
%
% Outputs: printed summary table, figures, challenge2_results.mat, and
% line-cut CSVs in the Challenge submission format.
%
% No toolboxes required - plain MATLAB on any CPU.

close all; clc; clear;

%% ===================== Configuration =====================
quick_mode = false;          % true = subsets [19 41], one strain window
use_batch  = true;           % batched IC-GN (dic_engine_batch): same math,
                             % all CPU cores, minutes instead of hours.
                             % false = force the sequential dic_engine.

if quick_mode
    subset_sweep    = [19 41];
    strain_windows  = 50;
else
    subset_sweep    = [9 19 29 39 49 59];   % px (odd); paper used 9..59 step 10
    strain_windows  = [20 50 80];           % PLS strain window sweep, px
end

band_half_disp   = 100;      % displacement: fit rows = center +- this many px
band_half_strain = 50;       % strain: fit rows with a FULL strain window

% engine parameters shared by every run (subset size is set per run)
P.grid_step    = 1;          % REQUIRED: 1 px line-cut spacing (no resampling)
P.search_range = 4;          % displacements are subpixel; seed search is cheap
P.black_thr    = 0.2;
P.zncc_min     = 0.5;
P.max_iter     = 30;
P.conv_tol     = 1e-3;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
challenge_dir = fullfile(this_dir, '..', '2D-Challenge 2.0');

% Star sets. lam0/lam1 define the linear period map lambda(x) over the
% image width (Reu et al. 2022, Eq. 1 / footnote 1). use_green = true for
% the strain sets (Lagrangian ground truth; for the first-harmonic
% amplitude the choice is immaterial - vy^2/2 only adds DC and 2nd
% harmonic). Stars 1-4 are disabled by default: their period map is
% OCR-ambiguous in our source (10->150 per footnote; VERIFY against the
% original PDF before enabling) and the paper's published clusters are
% for Stars 5-6 anyway.
stars = struct( ...
    'name',    {'Star1', 'Star2', 'Star3', 'Star4', 'Star5', 'Star6'}, ...
    'enabled', {false,   false,   false,   false,   true,    true}, ...
    'qoi',     {'disp',  'disp',  'strain','strain','disp',  'strain'}, ...
    'amp',     {0.5,     0.5,     0.05,    0.05,    0.5,     0.05}, ...
    'lam0',    {10,      10,      10,      10,      10,      10}, ...
    'lam1',    {150,     150,     150,     150,     300,     300}, ...
    'folder',  {fullfile('Star1NoNoise', 'Star1NoNoise'), ...
                fullfile('Star2Noise', 'Star2Noise'), ...
                fullfile('Star3NoNoiseStrain', 'Star3NoNoiseStrain'), ...
                fullfile('Star4NoiseStrain', 'Star4NoiseStrain'), ...
                fullfile('Star5LargeNoisy', 'Star5LargeNoisy'), ...
                fullfile('Star6StrainNoisy', 'Star6StrainNoisy')}, ...
    'ref',     {'DIC_Challenge_Wave_Reference_NoiseFree.tif', ...
                'DIC_Challenge_Wave_Reference_Noisy.tif', ...
                'DIC_Challenge_Wave_Reference_Strain.tif', ...
                'DIC_Challenge_Wave_Reference_Strain_Noise.tif', ...
                'DIC_Challenge_Star_Noise_Ref.tif', ...
                'DIC_Challenge_Star_Strain_Noise_Ref.tif'}, ...
    'floor_img', {'', ...                                    % Star 1: no noise-floor image
                'DIC_Challenge_Wave_Noisy2.tif', ...
                '', ...                                      % Star 3: no noise-floor image
                'DIC_Challenge_Wave_Strain_Noise2.tif', ...
                'DIC_Challenge_Star_Noise.tif', ...
                'DIC_Challenge_Star_Strain_Noise.tif'}, ...
    'def',     {'DIC_Challenge_Wave_Deformed_NoiseFree.tif', ...
                'DIC_Challenge_Wave_Deformed_Noisy.tif', ...
                'DIC_Challenge_Wave_Deformed_Strain.tif', ...
                'DIC_Challenge_Wave_Deformed_Strain_Noise.tif', ...
                'DIC_Challenge_Star_Noise_Def.tif', ...
                'DIC_Challenge_Star_Strain_Noise_Def.tif'});

% published per-code results (optional overlay): the two Results folders
% shipped with the Challenge download. NOTE the folder names are SWAPPED
% relative to their content (verified against the constant-amplitude
% ground truth): "DisplacementResults_Star6" holds the Star 5 displacement
% amplitudes, "StrainResults_Star5" holds the Star 6 strain amplitudes.
published = struct( ...
    'disp',   fullfile(challenge_dir, 'Results', 'DisplacementResults_Star6', 'Local.Affine.B.csv'), ...
    'strain', fullfile(challenge_dir, 'Results', 'StrainResults_Star5',       'Local.Affine.B.csv'));

%% ===================== Preflight =====================
addpath(fullfile(this_dir, '..', '6_DIC'));   % dic_engine(_batch).m, pls_strain.m

for s = find([stars.enabled])
    for f = {stars(s).ref, stars(s).def}
        assert(isfile(fullfile(challenge_dir, stars(s).folder, f{1})), ...
            'Missing %s image: %s (download from sem.org/dicchallenge, see README.md)', ...
            stars(s).name, f{1});
    end
end

if use_batch
    engine = @dic_engine_batch;
    % equivalence check on the REAL images: both engines solve the same
    % IC-GN problem; on a decimated grid of the first enabled star the
    % results must agree to within the convergence tolerance scale.
    s1 = find([stars.enabled], 1);
    assert(~isempty(s1), 'No star sets enabled.');
    refc = read_gray(fullfile(challenge_dir, stars(s1).folder, stars(s1).ref));
    defc = read_gray(fullfile(challenge_dir, stars(s1).folder, stars(s1).def));
    Pc = P; Pc.subset_size = 21; Pc.grid_step = 200;
    hbc  = floor(Pc.subset_size/2) + Pc.search_range + 2;
    icc  = round((size(refc,1) + 1) / 2);
    band = (icc - hbc):(icc + hbc);
    Rs = dic_engine(refc(band, :), defc(band, :), Pc);        % sequential
    Rb = dic_engine_batch(refc(band, :), defc(band, :), Pc);  % batched
    mm = isfinite(Rs.U) & isfinite(Rb.U);
    assert(nnz(mm) >= 10, 'equivalence check: too few common points');
    dmax = max(max(abs(Rs.U(mm) - Rb.U(mm))), max(abs(Rs.V(mm) - Rb.V(mm))));
    fprintf('Batch/sequential equivalence on %s: %d points, max |diff| = %.2e px\n', ...
        stars(s1).name, nnz(mm), dmax);
    assert(dmax < 0.02, 'batch engine disagrees with dic_engine (%.3g px)', dmax);
    if dmax > 5e-3
        warning('batch/sequential difference %.2g px is above the expected conv_tol scale.', dmax);
    end
else
    engine = @dic_engine;
end

%% ===================== Benchmark loop =====================
res  = struct('star', {}, 'qoi', {}, 'subset', {}, 'win_px', {}, ...
              'n', {}, 'ell10', {}, 'mei', {}, 'nfit', {});
cuts = struct('star', {}, 'qoi', {}, 'subset', {}, 'win_px', {}, ...
              'x', {}, 'q_floor', {}, 'q_def', {});
att  = struct('star', {}, 'qoi', {}, 'subset', {}, 'win_px', {}, ...
              'lam', {}, 'ratio', {}, 'l10', {});
t_all = tic;

for s = find([stars.enabled])
    st = stars(s);
    is_strain = strcmp(st.qoi, 'strain');
    use_green = is_strain;                    % Lagrangian ground truth

    ref = read_gray(fullfile(challenge_dir, st.folder, st.ref));
    def = read_gray(fullfile(challenge_dir, st.folder, st.def));
    flo = [];
    if ~isempty(st.floor_img)
        flo = read_gray(fullfile(challenge_dir, st.folder, st.floor_img));
    end
    [Himg, Wimg] = size(ref);
    lam = @(x) st.lam0 + (st.lam1 - st.lam0) / Wimg * (x - 1);
    ic = round((Himg + 1) / 2);               % center row of the image

    fprintf('\n=== %s (%s, %dx%d) ===\n', st.name, st.qoi, Himg, Wimg);

    for subset = subset_sweep
        P.subset_size = subset;
        half   = floor(subset / 2);
        margin = half + P.search_range + 2;
        if is_strain
            hb = margin + ceil(max(strain_windows)/2) + band_half_strain;
        else
            hb = margin + band_half_disp;
        end
        rows = (ic - hb):(ic + hb);
        assert(rows(1) >= 1 && rows(end) <= Himg, 'band exceeds image');

        t0 = tic;
        Rd = engine(ref(rows, :), def(rows, :), P);
        Rf = [];
        if ~isempty(flo)
            Rf = engine(ref(rows, :), flo(rows, :), P);
        end
        icb = find(Rd.gy == hb + 1, 1);       % grid row at the image center line
        assert(~isempty(icb), 'center row not on the DIC grid');
        lamx = lam(Rd.gx);                    % period at every grid column

        if ~is_strain
            % ---- displacement: noise from the center-row line-cut,
            %      attenuation from the per-column harmonic amplitude ----
            v_def = Rd.V(icb, :);
            v_flo = nan(size(v_def));
            if ~isempty(Rf), v_flo = Rf.V(icb, :); end
            n = std_finite(v_flo);
            A = harmonic_amp(Rd.V, Rd.gy, lamx);
            r = A / st.amp;
            l10 = ell10_fit(lamx, r);
            mei = n * l10;
            res(end+1) = struct('star', st.name, 'qoi', 'disp', 'subset', subset, ...
                'win_px', NaN, 'n', n, 'ell10', l10, 'mei', mei, 'nfit', nnz(isfinite(r))); %#ok<SAGROW>
            cuts(end+1) = struct('star', st.name, 'qoi', 'disp', 'subset', subset, ...
                'win_px', NaN, 'x', Rd.gx, 'q_floor', v_flo, 'q_def', v_def); %#ok<SAGROW>
            att(end+1) = struct('star', st.name, 'qoi', 'disp', 'subset', subset, ...
                'win_px', NaN, 'lam', lamx, 'ratio', r, 'l10', l10); %#ok<SAGROW>
            fprintf('  subset %2d:            n = %8.5f px   ell10%% = %6.1f px   MEI = %8.4f   (%.0f s)\n', ...
                subset, n, l10, mei, toc(t0));
        else
            % ---- strain: eps_yy per strain window; the harmonic fit uses
            %      only rows whose strain window is complete ----
            fit_rows = abs(Rd.gy - (hb + 1)) <= band_half_strain;
            for win = strain_windows
                [~, eyy_d] = pls_pair(Rd, P.grid_step, win, use_green);
                e_def = eyy_d(icb, :);
                e_flo = nan(size(e_def));
                if ~isempty(Rf)
                    [~, eyy_f] = pls_pair(Rf, P.grid_step, win, use_green);
                    e_flo = eyy_f(icb, :);
                end
                n = std_finite(e_flo);
                A = harmonic_amp(eyy_d(fit_rows, :), Rd.gy(fit_rows), lamx);
                r = A / st.amp;
                l10 = ell10_fit(lamx, r);
                mei = n^2 * l10;
                % Sanity bound: a linear smoothing pipeline
                % cannot resolve better than its ideal transfer, i.e.
                % sinc for the subset (plane-fit VALUE = mean filter)
                % times the PLS slope transfer 3(sin u - u cos u)/u^3 for
                % the strain window. Measured ell10 below this bound
                % means bias energy at the signal frequency (large
                % displacement swings on the strain star) is masking the
                % true attenuation - treat that ell10 as optimistic.
                lb = linspace(30, 299.9, 20000);
                u1 = pi*subset ./ lb;  u2 = pi*win ./ lb;
                tb = (sin(u1)./u1) .* (3*(sin(u2) - u2.*cos(u2)) ./ u2.^3);
                bnd = lb(find(tb < 0.9, 1, 'last'));
                if isfinite(l10) && ~isempty(bnd) && l10 < bnd
                    fprintf(['      NOTE: ell10 = %.1f px is below the ideal-filter ' ...
                        'bound %.1f px (%+.1f%%) - bias-inflated amplitude at long ' ...
                        'periods; ell10 and MEI are optimistic.\n'], ...
                        l10, bnd, 100*(l10/bnd - 1));
                end
                res(end+1) = struct('star', st.name, 'qoi', 'strain', 'subset', subset, ...
                    'win_px', win, 'n', n, 'ell10', l10, 'mei', mei, 'nfit', nnz(isfinite(r))); %#ok<SAGROW>
                cuts(end+1) = struct('star', st.name, 'qoi', 'strain', 'subset', subset, ...
                    'win_px', win, 'x', Rd.gx, 'q_floor', e_flo, 'q_def', e_def); %#ok<SAGROW>
                att(end+1) = struct('star', st.name, 'qoi', 'strain', 'subset', subset, ...
                    'win_px', win, 'lam', lamx, 'ratio', r, 'l10', l10); %#ok<SAGROW>
                fprintf('  subset %2d, VSG %2d px: n = %8.2e      ell10%% = %6.1f px   MEI = %8.2e (%.0f s)\n', ...
                    subset, win, n, l10, mei, toc(t0));
            end
        end
    end
end
fprintf('\nTotal benchmark time: %.0f s\n', toc(t_all));

%% ===================== Summary (MEI = avg of 3 lowest, per paper) =====================
for qoi = {'disp', 'strain'}
    sel = res(strcmp({res.qoi}, qoi{1}));
    if isempty(sel), continue; end
    m = sort([sel.mei]);
    m = m(isfinite(m));
    if isempty(m), continue; end
    mei3 = mean(m(1:min(3, numel(m))));
    if strcmp(qoi{1}, 'disp')
        fprintf('Displacement MEI (avg of 3 lowest): %.4f px^2\n', mei3);
    else
        fprintf('Strain MEI (avg of 3 lowest): %.3e\n', mei3);
    end
end

%% ===================== Figures =====================
for qoi = {'disp', 'strain'}
    sel = find(strcmp({res.qoi}, qoi{1}));
    if isempty(sel), continue; end
    is_str = strcmp(qoi{1}, 'strain');
    rep_win = strain_windows(min(2, numel(strain_windows)));   % representative VSG

    % (1) attenuation curves: measured ratio vs period + the 10% threshold
    figure('Name', sprintf('Challenge 2.0: %s attenuation', qoi{1}), ...
           'Position', [60, 60, 1100, 600]);
    hold on;
    for i = sel
        a = att(i);
        if is_str && a.win_px ~= rep_win, continue; end
        plot(a.lam, a.ratio, 'DisplayName', sprintf('subset %d', a.subset));
    end
    yline(0.9, 'r:', 'DisplayName', '10% attenuation');
    grid on; legend('Location', 'southeast'); ylim([-0.1, 1.3]);
    xlabel('period \lambda, px'); ylabel('measured / ground-truth amplitude');
    title(sprintf('%s: attenuation vs period (raw ratio; ell10%% from 12th-order fit)', qoi{1}));

    % (2) noise vs spatial resolution (the paper's Fig. 9 / Fig. 12 axes)
    figure('Name', sprintf('Challenge 2.0: %s noise vs resolution', qoi{1}));
    loglog([res(sel).ell10], [res(sel).n], 'o', 'LineWidth', 1.4, ...
        'DisplayName', 'LiteDIC');
    hold on; grid on;
    pub = published_points(published.(qoi{1}), ternary(is_str, 0.05, 0.5));
    if ~isempty(pub)
        loglog(pub(:,2), pub(:,1), 's--', 'LineWidth', 1.1, ...
            'DisplayName', 'Local.Affine.B (published)');
    end
    xlabel('spatial resolution l_{10%}, px');   % (\ell is not MATLAB TeX)
    ylabel(ternary(is_str, 'strain noise n (1\sigma)', 'displacement noise n (1\sigma), px'));
    legend('Location', 'northeast');
    title(sprintf('%s: measurement resolution vs spatial resolution', qoi{1}));
end

%% ===================== Save =====================
save(fullfile(this_dir, 'challenge2_results.mat'), 'res', 'cuts', 'att', 'P', ...
    'subset_sweep', 'strain_windows');
export_linecuts(this_dir, cuts, 'disp');
export_linecuts(this_dir, cuts, 'strain');
fprintf('Saved challenge2_results.mat and line-cut CSVs to %s\n', this_dir);

%% ===================== Local functions =====================

function img = read_gray(fname)
% Read an image as grayscale double in [0..1] without toolboxes.
    raw = imread(fname);
    if ndims(raw) == 3 && size(raw,3) >= 3
        raw = 0.2989*double(raw(:,:,1)) + 0.5870*double(raw(:,:,2)) + 0.1140*double(raw(:,:,3));
        img = raw / 255;
        return;
    end
    if isinteger(raw)
        img = double(raw) / double(intmax(class(raw)));
    else
        img = double(raw);
        mx = max(img(:));
        if mx > 1, img = img / mx; end
    end
end

function [exx, eyy] = pls_pair(R, step, win_px, use_green)
% PLS strains for one engine result (thin wrapper, keeps the loop tidy).
    [exx, eyy] = pls_strain(R.U, R.V, step, win_px, use_green);
end

function n = std_finite(v)
% Standard deviation over finite entries (1 sigma noise floor).
    v = v(isfinite(v));
    if numel(v) < 10, n = NaN; else, n = std(v); end
end

function A = harmonic_amp(q, yrows, lam_x)
% First-harmonic amplitude of a vertically sinusoidal field, per column.
% q: nR x nX measured field over the band; yrows: the nR row coordinates
% (px); lam_x: 1 x nX vertical period per column. Fits q(:,k) with
% [sin(2*pi*y/lam), cos(2*pi*y/lam), 1] by least squares (phase-free, so
% the fan's x-dependent phase reference never enters) and returns
% hypot(a,b). NaN where fewer than 10 valid rows.
    [~, nX] = size(q);
    w  = isfinite(q);
    qz = q;  qz(~w) = 0;
    th = 2*pi * (yrows(:) ./ lam_x);             % nR x nX
    Sn = sin(th) .* w;
    Cs = cos(th) .* w;
    s11 = sum(Sn.*Sn,1); s12 = sum(Sn.*Cs,1); s13 = sum(Sn,1);
    s22 = sum(Cs.*Cs,1); s23 = sum(Cs,1);     s33 = sum(w,1);
    b1  = sum(Sn.*qz,1); b2  = sum(Cs.*qz,1); b3  = sum(qz,1);
    M3 = zeros(3,3,nX); rhs = zeros(3,1,nX);
    M3(1,1,:) = s11; M3(1,2,:) = s12; M3(1,3,:) = s13;
    M3(2,1,:) = s12; M3(2,2,:) = s22; M3(2,3,:) = s23;
    M3(3,1,:) = s13; M3(3,2,:) = s23; M3(3,3,:) = s33;
    rhs(1,1,:) = b1; rhs(2,1,:) = b2; rhs(3,1,:) = b3;
    A = nan(1, nX);
    ok = s33 >= 10;
    if any(ok)
        cf = pagemldivide(M3(:,:,ok), rhs(:,:,ok));
        A(ok) = hypot(reshape(cf(1,1,:), 1, []), reshape(cf(2,1,:), 1, []));
    end
    A(~isfinite(A)) = NaN;
end

function l10 = ell10_fit(lp, r)
% Spatial resolution ell_10% per the Challenge 2.0 method: 12th-order
% polynomial fit of the amplitude ratio vs period, then the largest
% period at which the fit first drops below 0.9. Validated against the
% published Local.Affine.B amplitude data (reproduces the paper's
% per-subset values).
    l10 = NaN;
    m = isfinite(r) & isfinite(lp);
    if nnz(m) < 200
        warning('ell10_fit: only %d valid columns - not enough for the fit.', nnz(m));
        return;
    end
    lp = lp(m); r = r(m);
    [c, ~, mu] = polyfit(lp, r, 12);             % centered/scaled fit
    lf = linspace(min(lp), max(lp), 20000);
    rf = polyval(c, lf, [], mu);
    idx = find(rf < 0.9);
    idx = idx(idx < numel(lf));                  % need a right neighbor
    if isempty(idx) || numel(idx) >= numel(lf) - 1
        warning('ell10_fit: no 10%% crossing inside the period range.');
        return;
    end
    i = max(idx);                                % largest period below 0.9
    l10 = lf(i) + (0.9 - rf(i)) * (lf(i+1) - lf(i)) / max(rf(i+1) - rf(i), 1e-12);
end

function pub = published_points(csvfile, amp)
% (n, ell10%) per subset from a published Challenge results CSV
% (amplitude-envelope format: "Subset,<sizes>,,<sizes>" header, then
% "Pixel,Noise...,,Deformed..." columns). Returns [] if unavailable.
    pub = [];
    if ~isfile(csvfile), return; end
    try
        M = readmatrix(csvfile, 'NumHeaderLines', 2);
        H = readcell(csvfile, 'Range', '1:1');
        x = M(:, 1);
        sizes = [];
        for k = 2:numel(H)
            val = H{k};
            if ischar(val) || isstring(val), val = str2double(val); end
            if isnumeric(val) && isscalar(val) && isfinite(val)
                sizes(end+1) = val; %#ok<AGROW>
            else
                break;   % empty separator column reached
            end
        end
        nsub = numel(sizes);
        lam = @(xx) 10 + (300 - 10) / 4000 * (xx - 1);
        for k = 1:nsub
            noise_col = M(:, 1 + k);
            amp_col   = M(:, 1 + nsub + 1 + k);   % skip separator column
            nn = noise_col(isfinite(noise_col));
            if numel(nn) < 10, continue; end
            l10 = ell10_fit(lam(x).', (amp_col / amp).');
            if ~isfinite(l10), continue; end
            pub(end+1, :) = [std(nn), l10]; %#ok<AGROW>
        end
    catch err
        warning('Could not parse published results %s: %s', csvfile, err.message);
        pub = [];
    end
end

function export_linecuts(folder, cuts, qoi)
% Write the line-cuts in the Challenge submission format: one CSV per QOI,
% header rows "Subset,<sizes>,,<sizes>" and "Pixel,Noise...,,Deformed...",
% NaN outside the measured range. For strain, one file per strain window.
    sel = cuts(strcmp({cuts.qoi}, qoi));
    if isempty(sel), return; end
    wins = unique([sel.win_px]);
    if all(isnan(wins)), wins = NaN; end
    for w = wins(:)'
        if isnan(w), grp = sel(isnan([sel.win_px]));
        else,        grp = sel([sel.win_px] == w);
        end
        if isempty(grp), continue; end
        xmax = max(cellfun(@max, {grp.x}));
        X = (1:xmax)';
        ncol = numel(grp);
        NZ = nan(xmax, ncol); DF = nan(xmax, ncol);
        for k = 1:ncol
            NZ(grp(k).x, k) = grp(k).q_floor;
            DF(grp(k).x, k) = grp(k).q_def;
        end
        if isnan(w), fname = fullfile(folder, sprintf('LiteDIC_linecuts_%s.csv', qoi));
        else,        fname = fullfile(folder, sprintf('LiteDIC_linecuts_%s_VSG%d.csv', qoi, w));
        end
        fid = fopen(fname, 'w');
        subs = sprintf('%d,', [grp.subset]);
        % empty separator column between the Noise and Deformed blocks,
        % matching the published Challenge results format
        fprintf(fid, 'Subset,%s,,%s\n', subs(1:end-1), subs(1:end-1));
        fprintf(fid, 'Pixel%s,%s\n', repmat(',Noise', 1, ncol), ...
            repmat(',Deformed', 1, ncol));
        fmt = ['%d', repmat(',%.6g', 1, ncol), ',', repmat(',%.6g', 1, ncol), '\n'];
        for i = 1:xmax
            fprintf(fid, fmt, X(i), NZ(i,:), DF(i,:));
        end
        fclose(fid);
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
