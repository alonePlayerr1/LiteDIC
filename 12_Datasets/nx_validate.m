function nx_validate(varargin)
% NX_VALIDATE  Validate the LiteDIC subpixel engine against a FINITE-ELEMENT
% ground truth: a plate with a hole in tension, solved in NX / Simcenter
% Nastran (SOL 101). This closes the roadmap item "FE comparison of the
% hole-plate strain fields" and gives a third, independent validation next to
% the DICe gold (9_Validation_DICe) and the DIC Challenge (10_Challenge2).
%
% Pipeline (no NX or camera needed at run time):
%   FE surface field (nx_plate_hole_field.csv, mm)
%     -> FIELD_FROM_NX  : map mm->px window, remove rigid mean, amplify
%     -> GEN_DATASET    : render an EXACT reference/deformed speckle pair whose
%                         dense ground truth IS the FE field
%     -> DIC_ENGINE     : measure the pair
%     -> compare        : bias / RMS of (measured - FE) displacement, error map
%
% The FE displacement is amplified (see field_from_nx) so a modest linear
% solve becomes a few-pixel test; the ground truth is amplified by the same
% factor, so the comparison stays exact.
%
%   nx_validate                      % defaults (20 px/mm, amplify 300)
%   nx_validate('Scale',25,'Amplify',500,'GridStep',10,'SubsetSize',41)
%
% Options (name/value):
%   'Scale'       px per mm                         (default 20)
%   'Amplify'     FE displacement multiplier        (default 300)
%   'WinMM'       [Wmm Hmm] imaged window, mm        (default [25 50])
%   'GridStep'    DIC grid step, px                  (default 10)
%   'SubsetSize'  DIC subset, px (odd)               (default 41)
%   'Noise'       sensor noise std (0 = clean)       (default 0.005)
%   'Seed'        speckle RNG seed                   (default 7)
%
% Output: console bias/RMS report, plus <12_Datasets>/nx_out/ with the image
% pair, the ground truth, and nx_validate_report.png (GT vs measured vs error).
% Base MATLAB only -- no toolboxes.

    O.Scale      = 20;
    O.Amplify    = 300;
    O.WinMM      = [25 50];
    O.GridStep   = 10;
    O.SubsetSize = 41;
    O.Noise      = 0.005;
    O.Seed       = 7;
    for a = 1:2:numel(varargin), O.(varargin{a}) = varargin{a+1}; end

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '6_DIC'));      % dic_engine, pls_strain
    csv    = fullfile(here, 'nx_plate_hole_field.csv');
    outdir = fullfile(here, 'nx_out');
    if ~isfile(csv)
        error('nx_validate:nocsv', 'FE field not found: %s', csv);
    end

    % ---- 1. FE field -> pixel-space handles ----
    F = field_from_nx(csv, struct('scale', O.Scale, 'win_mm', O.WinMM, ...
        'amplify', O.Amplify));

    % ---- 2. render the exact reference/deformed pair from the FE field ----
    gen_dataset(struct('outdir', outdir, 'n_pairs', 1, 'imsize', F.imsize, ...
        'field_fun', F, 'seed', O.Seed, 'noise_std', O.Noise, ...
        'dots_per_kpx', 10, 'bits', 16, 'csv', false));

    ref = read_gray(fullfile(outdir, 'pair_0001_ref.png'));
    cur = read_gray(fullfile(outdir, 'pair_0001_def.png'));
    G   = load(fullfile(outdir, 'pair_0001_gt.mat'));   % U,V,EXX,EYY,EXY

    % ---- 3. mask the hole: blacken the disc in both frames and exclude the
    %        hole + a subset-half rim (subsets straddling the hole cannot be
    %        measured -- exactly what a real DIC run does at a cut-out) ----
    [H, W]   = size(ref);
    [XX, YY] = meshgrid(1:W, 1:H);              % 1-based pixel centres
    hc = F.hole.center_px + 1;                  % 0-based -> 1-based
    dist = hypot(XX - hc(1), YY - hc(2));
    in_hole = dist <= F.hole.r_px;
    ref(in_hole) = 0;  cur(in_hole) = 0;        % a genuine black hole
    keep = dist > F.hole.r_px + O.SubsetSize/2; % rim subsets are unmeasurable

    % ---- 4. measure with the standard engine ----
    P.subset_size  = O.SubsetSize;
    P.grid_step    = O.GridStep;
    P.search_range = 15;
    P.black_thr    = 0.05;    % the blackened hole (0) is masked; bg (0.12) kept
    P.zncc_min     = 0.5;
    P.max_iter     = 30;
    P.conv_tol     = 1e-3;
    P.mask         = keep;    % drop hole + rim nodes (see above)
    R = dic_engine(ref, cur, P);

    % ---- 5. compare measured vs FE ground truth at the grid nodes ----
    % R.gx / R.gy are 1-based pixel indices; the dense GT element (r,c) is the
    % displacement at pixel (row r, col c), so GT at a node is a direct index.
    Ugt = G.U(R.gy, R.gx);
    Vgt = G.V(R.gy, R.gx);
    m   = ~isnan(R.U) & ~isnan(R.V) & isfinite(Ugt) & isfinite(Vgt);
    du  = R.U(m) - Ugt(m);
    dv  = R.V(m) - Vgt(m);

    bias_u = mean(du);           bias_v = mean(dv);
    rms_u  = sqrt(mean(du.^2));  rms_v  = sqrt(mean(dv.^2));

    BIAS_TOL = 0.02; RMS_TOL = 0.05;    % same criterion as the DICe check
    pass = max(abs([bias_u bias_v])) < BIAS_TOL && max([rms_u rms_v]) < RMS_TOL;

    fprintf('\n===== NX (FE) validation =====\n');
    fprintf('image %d x %d px @ %g px/mm, amplify %g, peak |v| %.2f px\n', ...
        F.imsize(2), F.imsize(1), O.Scale, O.Amplify, F.info.peak_v_px);
    fprintf('compared points: %d\n', nnz(m));
    fprintf('bias   u %.4f px   v %.4f px\n', bias_u, bias_v);
    fprintf('RMS    u %.4f px   v %.4f px\n', rms_u, rms_v);
    fprintf('Criterion |bias| < %.2f px and RMS < %.2f px  ->  %s\n', ...
        BIAS_TOL, RMS_TOL, ternary(pass, 'PASS', 'FAIL'));

    % ---- 6. report figure: GT v, measured v, error map ----
    err = nan(size(R.V)); err(m) = dv;
    fig = figure('Visible', 'off', 'Position', [100 100 1100 420], 'Color', 'w');
    subplot(1,3,1); imagesc(R.gx, R.gy, Vgt); axis image; colorbar;
        title('FE v (ground truth), px'); set(gca,'YDir','reverse');
    subplot(1,3,2); imagesc(R.gx, R.gy, R.V); axis image; colorbar;
        title('LiteDIC v (measured), px'); set(gca,'YDir','reverse');
    subplot(1,3,3); imagesc(R.gx, R.gy, err); axis image; colorbar;
        title(sprintf('v error (meas-FE), RMS %.3f px', rms_v));
        set(gca,'YDir','reverse');
    print(fig, fullfile(outdir, 'nx_validate_report.png'), '-dpng', '-r120');
    close(fig);
    fprintf('report -> %s\n', fullfile(outdir, 'nx_validate_report.png'));
end

% ---- helpers (toolbox-free) ----
function g = read_gray(pathname)
    I = imread(pathname);
    if isinteger(I), s = double(intmax(class(I))); else, s = 1; end
    if ndims(I) == 3 && size(I,3) >= 3
        I = double(I);
        g = (0.2989*I(:,:,1) + 0.5870*I(:,:,2) + 0.1140*I(:,:,3)) / s;
    else
        g = double(I) / s;
    end
    mx = max(g(:)); if mx > 1, g = g / mx; end
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
