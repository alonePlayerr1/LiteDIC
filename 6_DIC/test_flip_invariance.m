function out = test_flip_invariance(varargin)
% TEST_FLIP_INVARIANCE  Prove the 2D DIC engine is orientation-agnostic and
% expose that the *presentation* layer is not.
%
%   out = test_flip_invariance()
%   out = test_flip_invariance('Frames',{'plate_def_02.tiff'})
%   out = test_flip_invariance('Frames',{'plate_def_11.tiff'},'Search',25)
%
% WHY THIS TEST EXISTS
%   A correct DIC pipeline must be EQUIVARIANT under a rigid image
%   transform: if the same transform T is applied to the reference and the
%   deformed frame, the recovered field is the same physical field, merely
%   re-expressed in the transformed coordinate frame. Consequently every
%   SCALAR INVARIANT of that field - |displacement|, tr(e)=exx+eyy, the
%   principal strains {e1,e2}, von Mises and max shear - must be IDENTICAL
%   (up to interpolation noise) no matter how the images are flipped or
%   rotated. The COMPONENTS transform as a vector / tensor (under 90 deg
%   exx<->eyy swap; under a mirror exy flips sign).
%
%   The test applies the 8 dihedral (D4) transforms - pixel-exact flips and
%   90 deg rotations, NO resampling - and splits them into two classes:
%     SWAP class  (transpose, antitranspose): pure axis relabel, reverses
%                 no direction. This isolates "does the engine have a baked-
%                 in up?" -> it should reproduce the baseline to the
%                 interpolation floor.
%     REVERSAL    (rot90/180/270, fliplr, flipud): reverse one or both
%                 axes. Any residual here is directional SUBPIXEL BIAS (the
%                 S-curve, cf. bias_scurve.m), which is real, expected and
%                 grows with the displacement gradient - so it should
%                 shrink toward the swap floor on a low-load frame.
%
%   Sweeping frames from low to high load therefore separates the two
%   questions cleanly: SWAP ~ floor at every load = engine is orientation-
%   free; REVERSAL growing with load = benign subpixel bias, not a bug.
%
%   The run also renders the eyy panel per orientation exactly as
%   run_series_2d.m / make_eyy_figure.m plot it. Those panels DIFFER for the
%   same physical state - that is the presentation bug: the plotting layer
%   hardcodes "the strain = eyy" instead of reading the load axis from data.
%
% Base MATLAB only - no toolboxes, no GPU. Uses dic_engine.m + pls_strain.m.

%% -------- options --------
here = fileparts(mfilename('fullpath'));
O.Ref       = fullfile(here, 'plate_ref.tiff');
O.Frames    = {'plate_def_02.tiff','plate_def_06.tiff','plate_def_11.tiff'};
O.GridStep  = 20;
O.SubsetSize= 41;
O.Search    = 15;
O.OutDir    = fullfile(here, 'flip_test_out');
for a = 1:2:numel(varargin), O.(varargin{a}) = varargin{a+1}; end
if ischar(O.Frames), O.Frames = {O.Frames}; end
if ~exist(O.OutDir, 'dir'), mkdir(O.OutDir); end

P.subset_size  = O.SubsetSize;
P.grid_step    = O.GridStep;
P.search_range = O.Search;
P.black_thr    = 0.2;
P.zncc_min     = 0.5;
P.max_iter     = 30;
P.conv_tol     = 1e-3;
strain_win = max(50, 3 * P.grid_step);

ref0 = load_gray(O.Ref);
[H, W] = size(ref0);
fprintf('Flip-invariance sweep\n  ref : %s  (%d rows x %d cols)\n\n', O.Ref, H, W);

%% -------- sweep frames --------
% SWAP class = pure axis relabel; REVERSAL class = axis reversal(s).
swapNames = {'transpose','antitranspose'};
hdr = sprintf('%-16s %9s %9s %11s %11s %8s', 'frame', 'max|d|px', 'RMSeyy', ...
    'SWAPfloor', 'REVworst', 'REV/base');
fprintf('%s\n', hdr); fprintf('%s\n', repmat('-', 1, numel(hdr)));

sweep = struct('frame', {}, 'maxdisp', {}, 'rmseyy', {}, 'swapfloor', {}, ...
    'revworst', {}, 'revrel', {}, 'res', {}, 'panels', {}, 'Ib', {});
for q = 1:numel(O.Frames)
    def0 = load_gray(fullfile(here, O.Frames{q}));
    [res, panels, Ib] = run_one_pair(ref0, def0, P, strain_win, swapNames);

    baseScale = max(abs(Ib.dispmag(:)));            % field magnitude (px)
    rmseyy    = sqrt(mean(Ib.trE(:).^2, 'omitnan')); % rough strain scale proxy
    isSwap = ismember({res.name}, swapNames);
    dsw = arrayfun(@(r) r.metrics.dispmag(2), res(isSwap));
    drv = arrayfun(@(r) r.metrics.dispmag(2), res(~isSwap & ~strcmp({res.name},'identity')));
    swapfloor = max(dsw);
    revworst  = max(drv);
    revrel    = revworst / max(baseScale, eps);

    fprintf('%-16s %9.3f %9.2e %11.2e %11.2e %7.1f%%\n', O.Frames{q}, ...
        baseScale, rmseyy, swapfloor, revworst, 100*revrel);
    sweep(q) = struct('frame', O.Frames{q}, 'maxdisp', baseScale, ...
        'rmseyy', rmseyy, 'swapfloor', swapfloor, 'revworst', revworst, ...
        'revrel', revrel, 'res', res, 'panels', {panels}, 'Ib', Ib);
end

%% -------- verdict --------
% Engine is orientation-free iff the SWAP class sits at the interpolation
% floor at EVERY load (independent of displacement magnitude). The REVERSAL
% class is graded separately as a quality (bias) metric, not a pass gate.
swapMax = max([sweep.swapfloor]);
fprintf('\n--- interpretation ---\n');
fprintf('SWAP class (pure axis relabel) worst over all frames: %.2e px\n', swapMax);
if swapMax < 2e-2
    fprintf(['ENGINE: orientation-agnostic. A pure axis swap (the horizontal-\n' ...
        '        vs-vertical mounting case) reproduces the field to the\n' ...
        '        interpolation floor at every load.\n']);
else
    fprintf(['ENGINE: CHECK - pure axis swap should be at the interp floor but\n' ...
        '        is not; this would indicate a real handedness asymmetry.\n']);
end
% does the reversal-class bias grow with load (=> subpixel bias, benign)?
if numel(sweep) >= 2
    lo = sweep(1); hi = sweep(end);
    fprintf(['REVERSAL bias vs load: %.2e px @ max|d|=%.2f  ->  %.2e px @ max|d|=%.2f\n' ...
        '        %s\n'], lo.revworst, lo.maxdisp, hi.revworst, hi.maxdisp, ...
        tern(hi.revworst > 2*lo.revworst, ...
        'grows with load => directional subpixel bias (S-curve), expected/benign.', ...
        'roughly flat with load => dominated by the interpolation floor.'));
end

%% -------- detailed table + figure for the highest-load frame --------
last = sweep(end);
fprintf('\nDetailed per-transform RMS (frame %s):\n', last.frame);
fprintf('%-14s %10s %11s %11s %11s %11s %11s\n', 'transform', 'validpts', ...
    '|disp|px', 'tr(e)', 'e1', 'e2', 'vonMises');
for t = 1:numel(last.res)
    m = last.res(t).metrics;
    fprintf('%-14s %10d %11.2e %11.2e %11.2e %11.2e %11.2e\n', last.res(t).name, ...
        last.res(t).valid, m.dispmag(2), m.trE(2), m.e1(2), m.e2(2), m.vm(2));
end

show = {'identity','rot90','flipud','transpose'};
f = figure('Visible', 'off', 'Position', [50 50 1400 380]);
for s = 1:numel(show)
    t = find(strcmp({last.res.name}, show{s}), 1);
    subplot(1, numel(show), s);
    h = imagesc(last.panels{t}.gx, last.panels{t}.gy, last.panels{t}.eyy);
    set(h, 'AlphaData', isfinite(last.panels{t}.eyy));
    axis image; colorbar;
    title(sprintf('\\epsilon_{yy} panel  [%s]', show{s}));
end
sgtitle(['Same physical strain state, one hardcoded panel (\epsilon_{yy}): ' ...
    'the plot changes with orientation - presentation bug']);
print(f, fullfile(O.OutDir, 'presentation_bug_eyy.png'), '-dpng', '-r120');
close(f);
fprintf('\nSaved figure: %s\n', fullfile(O.OutDir, 'presentation_bug_eyy.png'));

out.sweep = sweep;  out.P = P;
end

% ===================== per-pair (8 D4 transforms) =====================

function [res, panels, Ib] = run_one_pair(ref0, def0, P, strain_win, ~)
[H, W] = size(ref0);
T = { ...
    'identity',    @(A) A; ...
    'rot90',       @(A) rot90(A, 1); ...
    'rot180',      @(A) rot90(A, 2); ...
    'rot270',      @(A) rot90(A, 3); ...
    'fliplr',      @(A) fliplr(A); ...
    'flipud',      @(A) flipud(A); ...
    'transpose',   @(A) A.'; ...
    'antitranspose', @(A) rot90(fliplr(A), 1) };
nT = size(T, 1);
[Col0, Row0] = meshgrid(1:W, 1:H);

Rb = dic_engine(ref0, def0, P);
[bx, by] = meshgrid(Rb.gx, Rb.gy);
Ib = invariants(Rb, P.grid_step, strain_win);
fld = fieldnames(Ib);

res = struct('name', {}, 'valid', {}, 'metrics', {});
panels = cell(nT, 1);
for t = 1:nT
    op = T{t,2};
    reft = op(ref0);  deft = op(def0);
    RowT = op(Row0);  ColT = op(Col0);

    R = dic_engine(reft, deft, P);
    [I, comp] = invariants(R, P.grid_step, strain_win);
    panels{t} = struct('gx', R.gx, 'gy', R.gy, 'eyy', comp.eyy);

    [gxt, gyt] = meshgrid(R.gx, R.gy);
    lin = sub2ind(size(reft), round(gyt(:)), round(gxt(:)));
    oc = ColT(lin);  orr = RowT(lin);

    m = struct();
    for f = 1:numel(fld)
        key = fld{f};  vals = I.(key);
        good = isfinite(vals(:)) & isfinite(oc) & isfinite(orr);
        if nnz(good) < 10, m.(key) = [NaN NaN]; continue; end
        F = scatteredInterpolant(oc(good), orr(good), vals(good), 'linear', 'none');
        d = F(bx, by) - Ib.(key);
        cm = isfinite(d) & isfinite(Ib.(key));
        m.(key) = [max(abs(d(cm))), sqrt(mean(d(cm).^2))];
    end
    res(t) = struct('name', T{t,1}, 'valid', nnz(isfinite(comp.eyy)), 'metrics', m);
end
end

% ===================== helpers =====================

function [I, comp] = invariants(R, step, win_px)
    [exx, eyy, exy] = pls_strain(R.U, R.V, step, win_px, false);
    comp.exx = exx;  comp.eyy = eyy;  comp.exy = exy;
    I.dispmag = sqrt(R.U.^2 + R.V.^2);
    I.trE     = exx + eyy;
    a = (exx + eyy) / 2;
    b = sqrt(((exx - eyy) / 2).^2 + exy.^2);
    I.e1 = a + b;  I.e2 = a - b;
    I.vm = sqrt(I.e1.^2 - I.e1.*I.e2 + I.e2.^2);
    I.maxshear = b;
end

function g = load_gray(fn)
    A = imread(fn);
    if ndims(A) == 3, A = mean(double(A), 3); else, A = double(A); end
    g = A;
    if max(g(:)) > 1.5, g = g / max(g(:)); end
end

function s = tern(c, a, b), if c, s = a; else, s = b; end, end
