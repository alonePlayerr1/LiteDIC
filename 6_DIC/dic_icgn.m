% Subpixel 2D DIC engine: ZNSSD + IC-GN + RG-DIC + PLS strains. CPU only.
%
%   * ZNSSD correlation criterion - invariant to linear illumination
%     changes between frames;
%   * subpixel registration by the inverse-compositional Gauss-Newton
%     method (IC-GN, Baker & Matthews): the Hessian is assembled ONCE per
%     point, first-order (affine) subset shape function;
%   * spline interpolation of image intensity (registration accuracy
%     better than 0.01 px);
%   * reliability-guided initial guesses (RG-DIC): the solution grows from
%     a single seed point to its neighbors in order of correlation quality,
%     so each new point starts from a converged neighbor - 3-5 iterations
%     per point. The front is kept in a binary heap (O(N log N)), and a
%     point is finalized only on successful convergence, so a failure from
%     one neighbor does not block convergence from another;
%   * strains by pointwise least squares (PLS): displacements inside a
%     strain window are fitted with a plane and differentiated
%     ANALYTICALLY - noise is averaged before differentiation, unlike
%     finite differences. The window is specified in PIXELS, so the strain
%     field does not depend on the grid step. Strains are computed only at
%     valid nodes: extrapolating the LSQ plane into the hole would create
%     large fake strain peaks at the boundary;
%   * optional Green-Lagrange strain measure (for large rotations);
%   * built-in self-test: a known rigid subpixel shift must be recovered
%     to < 0.02 px with zero strains (guards against sign/axis errors that
%     are invisible in "nice-looking" fields);
%   * processes the full loading series against the reference frame; the
%     frames are discovered automatically by filename pattern.
%
% The numerical kernels live in dic_engine.m (correlation) and
% pls_strain.m (strains) in this folder, so other scripts (e.g. the
% Challenge 2.0 benchmark in 10_Challenge2) can call them on arbitrary
% inputs. This script is the acquisition + visualization front end.
%
% No toolboxes and no GPU required - plain MATLAB on any CPU.

close all; clc; clear;

%% ===================== Parameters =====================
P.subset_size  = 41;     % correlation subset, px (odd)
P.grid_step    = 10;     % grid step, px
P.search_range = 15;     % integer search ONLY for the seed point, px
P.black_thr    = 0.2;    % "black region" threshold (the hole)
P.zncc_min     = 0.5;    % minimum quality; points below are rejected
P.max_iter     = 30;     % IC-GN iteration cap
P.conv_tol     = 1e-3;   % convergence on ||dp||, px
P.strain_win_px = 50;    % PLS strain window, PIXELS (grid-step independent)
P.use_green    = false;  % true = Green-Lagrange strains (large rotations)

run_selftest = true;     % synthetic self-test before the real analysis

% Acquisition layer: the series is discovered by filename pattern, so any
% image set can be analyzed by changing these two lines (or the folder).
% ref_pattern must match exactly one file; def_pattern matches the deformed
% frames, processed in sorted (chronological) order.
ref_pattern = 'plate_ref.tiff';
def_pattern = 'plate_def_*.tiff';

%% ===================== Self-test (engine verification) =====================
if run_selftest
    fprintf('=== Self-test: rigid shift (+0.30, -0.20) px ===\n');
    rng(0);
    img = conv2(rand(220), gauss_kernel(2.0), 'same');
    img = img(11:210, 11:210);                       % crop edge effects
    img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    TU = 0.30; TV = -0.20;                           % true shift (x, y)
    Fi = griddedInterpolant({1:size(img,1), 1:size(img,2)}, img, 'spline', 'nearest');
    [XX, YY] = meshgrid(1:size(img,2), 1:size(img,1));
    shifted = reshape(Fi(YY(:) - TV, XX(:) - TU), size(img));
    R = dic_engine(img, shifted, P);
    [exx, eyy, exy] = pls_strain(R.U, R.V, P.grid_step, P.strain_win_px, P.use_green);
    eu = mean(R.U(:), 'omitnan') - TU;
    ev = mean(R.V(:), 'omitnan') - TV;
    smax = max([max(abs(exx(:))), max(abs(eyy(:))), max(abs(exy(:)))], [], 'omitnan');
    fprintf('  error u: %+.5f px, error v: %+.5f px, max|strain|: %.1e\n', eu, ev, smax);
    if abs(eu) < 0.02 && abs(ev) < 0.02 && smax < 5e-4
        fprintf('  SELF-TEST PASSED: no spurious strains.\n\n');
    else
        warning('SELF-TEST FAILED - check the parameters!');
    end
end

%% ===================== Acquisition layer =====================
this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
[ref_name, cur_names] = find_series(this_dir, ref_pattern, def_pattern);
ref_image = read_gray(fullfile(this_dir, ref_name));

nF = numel(cur_names);
series = struct('name', cur_names, 'mean_u', [], 'mean_v', [], 'zncc', [], ...
    'U', [], 'V', [], 'Q', [], 'exx', [], 'eyy', [], 'exy', []);

%% ============ Engine + strains over the whole series ============
for k = 1:nF
    cur_image = read_gray(fullfile(this_dir, cur_names{k}));
    if ~isequal(size(cur_image), size(ref_image))
        error('Size of %s does not match the reference frame.', cur_names{k});
    end
    t0 = tic;
    R = dic_engine(ref_image, cur_image, P);
    [exx, eyy, exy] = pls_strain(R.U, R.V, P.grid_step, P.strain_win_px, P.use_green);
    series(k).mean_u = mean(R.U(:), 'omitnan');
    series(k).mean_v = mean(R.V(:), 'omitnan');
    series(k).zncc   = mean(R.Q(:), 'omitnan');
    % full fields per frame - consumed by 9_Validation_DICe/validate_dice.m
    series(k).U = R.U; series(k).V = R.V; series(k).Q = R.Q;
    series(k).exx = exx; series(k).eyy = eyy; series(k).exy = exy;
    fprintf('%s -> %s: %d points, %.1f s\n', ref_name, cur_names{k}, ...
        nnz(~isnan(R.U)), toc(t0));
    fprintf('  mean u(x) = %+.3f px, v(y) = %+.3f px, ZNCC = %.4f, iterations %.2f\n', ...
        series(k).mean_u, series(k).mean_v, series(k).zncc, mean(R.NIT(:), 'omitnan'));
end
% R, exx... hold the LAST frame of the series - that is what gets plotted.

%% ===================== Visualization layer =====================
cur_image = read_gray(fullfile(this_dir, cur_names{end}));

figure('Name', 'DIC: input frames');
subplot(1,2,1); imagesc(ref_image); colormap(gray); axis image off;
title(sprintf('Reference frame: %s', ref_name), 'Interpreter', 'none');
subplot(1,2,2); imagesc(cur_image); colormap(gray); axis image off;
title(sprintf('Deformed: %s', cur_names{end}), 'Interpreter', 'none');

valid_mask = ~isnan(R.Q);
fields = {R.U, R.V, R.Q, exx, eyy, exy};
titles = {'Displacement u (X), px', 'Displacement v (Y), px', 'ZNCC quality', ...
          'Strain \epsilon_{xx} (PLS)', 'Strain \epsilon_{yy} (PLS)', ...
          'Strain \epsilon_{xy} (PLS)'};
figure('Name', 'DIC: fields (IC-GN + PLS)', 'Position', [80, 80, 1250, 800]);
for k = 1:6
    subplot(2,3,k);
    h = imagesc(R.gx, R.gy, fields{k});
    set(h, 'AlphaData', valid_mask & ~isnan(fields{k}));
    colorbar; axis image; title(titles{k}); xlabel('X, px'); ylabel('Y, px');
    if k >= 4   % strains: 1-99 percentile color scale (outliers don't crush it)
        lim = robust_lim(fields{k});
        if lim(2) > lim(1), clim(lim); end
    end
end

figure('Name', 'DIC: vector field');
imagesc(cur_image); colormap(gray); axis image; hold on;
[Xg, Yg] = meshgrid(R.gx, R.gy);
vi = find(valid_mask);
quiver(Xg(vi), Yg(vi), R.U(vi), R.V(vi), 2, 'r', 'LineWidth', 1.1);
title('Displacement vector field (last frame)'); xlabel('X'); ylabel('Y'); hold off;

figure('Name', 'DIC: evolution over the loading series');
frames = 1:nF;
yyaxis left;
plot(frames, abs([series.mean_v]), 'o-', 'LineWidth', 1.5);
ylabel('|mean v|, px');
yyaxis right;
plot(frames, [series.zncc], 's--', 'LineWidth', 1.5);
ylabel('mean ZNCC'); ylim([0.9, 1.0]);
set(gca, 'XTick', frames, 'XTickLabel', strrep(cur_names, '_', '\_'));
xlabel('Series frame'); grid on;
title('Displacement growth over the loading series');

save(fullfile(this_dir, 'dic_results.mat'), 'R', 'exx', 'eyy', 'exy', 'series', 'P');
fprintf('\nResults saved to dic_results.mat\n');

%% ===================== Local functions =====================
% (the numerical kernels dic_engine.m and pls_strain.m live in their own
%  files in this folder; only script-specific helpers remain here)

function lim = robust_lim(f)
% Color limits at the 1st and 99th percentiles (no toolboxes needed).
    v = sort(f(isfinite(f)));
    if isempty(v), lim = [0 1]; return; end
    n = numel(v);
    lim = [v(max(1, round(0.01*n))), v(min(n, round(0.99*n)))];
end

function K = gauss_kernel(sigma)
% Gaussian kernel for the synthetic speckle in the self-test.
    r = ceil(3*sigma);
    [x, y] = meshgrid(-r:r);
    K = exp(-(x.^2 + y.^2) / (2*sigma^2));
    K = K / sum(K(:));
end

function [ref_name, cur_names] = find_series(folder, ref_pattern, def_pattern)
% Discover the image series in a folder: one reference frame plus the
% deformed frames in sorted (chronological) order. Fails with a clear
% message instead of erroring later inside imread.
    d = dir(fullfile(folder, ref_pattern));
    if isempty(d)
        error('Reference frame not found: %s in %s', ref_pattern, folder);
    elseif numel(d) > 1
        error('Reference pattern %s matches %d files in %s - it must match exactly one.', ...
            ref_pattern, numel(d), folder);
    end
    ref_name = d(1).name;
    d = dir(fullfile(folder, def_pattern));
    cur_names = setdiff(sort_frames({d.name}), {ref_name}, 'stable');
    if isempty(cur_names)
        error('No deformed frames found: %s in %s', def_pattern, folder);
    end
    fprintf('Series: %s + %d deformed frames (%s ... %s)\n', ...
        ref_name, numel(cur_names), cur_names{1}, cur_names{end});
end

function g = to_gray(img)
% Grayscale conversion without the Image Processing Toolbox.
    if ndims(img) == 3 && size(img,3) >= 3
        img = double(img);
        g = 0.2989*img(:,:,1) + 0.5870*img(:,:,2) + 0.1140*img(:,:,3);
    else
        g = double(img);
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

function names = sort_frames(names)
% Chronological sort of frame names. Real cameras/DAQs often number frames
% WITHOUT zero-padding (img_2.jpg ... img_10.jpg): a plain lexicographic
% sort puts img_10 before img_2 and the series is processed out of order.
% When every name carries a trailing number before the extension, sort by
% that number; otherwise fall back to the lexicographic sort. Zero-padded
% series sort identically either way, so existing datasets are unaffected.
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
