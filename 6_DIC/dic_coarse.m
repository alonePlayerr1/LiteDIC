% 2D DIC (digital image correlation) - fast integer-pixel version, CPU only.
%
% Speed idea: for every integer shift (du,dv) the NCC field is computed for
% the WHOLE grid at once with box filters (separable convolutions) instead
% of nested per-point loops - orders of magnitude faster than the naive
% element-wise search. Subpixel refinement: parabola fit around the NCC peak.
%
% This is the quick-look tool (mean NCC ~0.95 on the test pair). For
% subpixel-accurate fields use dic_icgn.m in the same folder.
close all; clc; clear;

% Parameters
subset_size  = 41;       % subset size, px (odd)
grid_step    = 10;       % grid step, px (5 gives a denser grid, a bit slower)
search_range = 15;       % search range, px
black_threshold = 0.2;   % "black region" threshold (the hole)

% Load images (from the script's own folder). The pair is discovered by
% pattern: the reference frame plus the LAST deformed frame of the series
% (largest deformation = the interesting quick look).
ref_pattern = 'plate_ref.tiff';
def_pattern = 'plate_def_*.tiff';

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
d_ref = dir(fullfile(this_dir, ref_pattern));
d_def = dir(fullfile(this_dir, def_pattern));
def_names = setdiff(sort_frames({d_def.name}), {d_ref.name}, 'stable');
if isempty(d_ref) || isempty(def_names)
    error('Image series not found: %s / %s in %s.', ref_pattern, def_pattern, this_dir);
end
fprintf('Pair: %s -> %s\n', d_ref(1).name, def_names{end});
ref_image = read_gray(fullfile(this_dir, d_ref(1).name));
def_image = read_gray(fullfile(this_dir, def_names{end}));
image_size = size(ref_image);
if ~isequal(size(def_image), image_size)
    error('Image sizes do not match.');
end

% Show the input images
figure;
subplot(1,2,1); imagesc(ref_image); colormap(gray); axis image; axis off;
title('Reference image');
subplot(1,2,2); imagesc(def_image); colormap(gray); axis image; axis off;
title('Deformed image');

half = floor(subset_size/2);
pad  = search_range + half;
N    = subset_size^2;
Rp = sym_pad(ref_image, pad);
Dp = sym_pad(def_image, pad);

% Mask of "black" subsets (max brightness within the window)
Rmax = movmax(movmax(Rp, subset_size, 1), subset_size, 2);

% Vectorized CPU computation
fprintf('Vectorized CPU computation.\n');
kcol = ones(subset_size, 1); krow = ones(1, subset_size);
Rg = Rp; Dg = Dp;

boxsum = @(X) conv2(conv2(X, kcol, 'same'), krow, 'same');   % window sum

Rmean = boxsum(Rg) / N;
Rstd  = sqrt(max(boxsum(Rg.^2)/N - Rmean.^2, 0));

% Analysis grid (in padded-image coordinates)
gy = (pad+half+1) : grid_step : (pad+image_size(1)-half);
gx = (pad+half+1) : grid_step : (pad+image_size(2)-half);
nGy = numel(gy); nGx = numel(gx);
S = 2*search_range + 1;

% NCC cube: per grid point, the NCC map over all shifts
tic;
cube = -2 * ones(nGy, nGx, S, S);
for iu = 1:S
    du = iu - search_range - 1;       % row shift (Y)
    for iv = 1:S
        dv = iv - search_range - 1;   % column shift (X)
        Dsh = circshift(Dg, [-du, -dv]);
        Dmean = boxsum(Dsh) / N;
        Dstd  = sqrt(max(boxsum(Dsh.^2)/N - Dmean.^2, 0));
        cross = boxsum(Rg.*Dsh)/N - Rmean.*Dmean;
        ncc = cross ./ (Rstd .* Dstd);
        ncc(~isfinite(ncc)) = -2;
        cube(:, :, iu, iv) = ncc(gy, gx);
    end
end
fprintf('Correlation: %.2f s  (%d shifts, grid %dx%d)\n', toc, S*S, nGy, nGx);

% Integer peak per grid point
flat = reshape(cube, nGy, nGx, S*S);
[ncc_values, idx] = max(flat, [], 3);
[pu, pv] = ind2sub([S, S], idx);
du0 = pu - search_range - 1;   % vertical (Y)
dv0 = pv - search_range - 1;   % horizontal (X)

% Subpixel correction (parabola around the peak) - vectorized
[I, J] = ndgrid(1:nGy, 1:nGx);
val = @(a,b) cube(sub2ind([nGy,nGx,S,S], I, J, min(max(a,1),S), min(max(b,1),S)));
c0 = val(pu, pv);
sub_du = zeros(nGy, nGx); sub_dv = zeros(nGy, nGx);
inU = pu>1 & pu<S; inV = pv>1 & pv<S;
fm = val(pu-1,pv); fp = val(pu+1,pv); d = fm - 2*c0 + fp;
t = 0.5*(fm - fp)./d; t(d==0 | ~isfinite(t) | abs(t)>1) = 0; sub_du(inU) = t(inU);
fm = val(pu,pv-1); fp = val(pu,pv+1); d = fm - 2*c0 + fp;
t = 0.5*(fm - fp)./d; t(d==0 | ~isfinite(t) | abs(t)>1) = 0; sub_dv(inV) = t(inV);

displacement_y = du0 + sub_du;   % vertical (rows, Y)
displacement_x = dv0 + sub_dv;   % horizontal (columns, X)

% Validity mask
refmax_grid = Rmax(gy, gx);
valid = (refmax_grid >= black_threshold) & (ncc_values > -1);
displacement_x(~valid) = NaN;
displacement_y(~valid) = NaN;
ncc_values(~valid) = NaN;

% Strains (finite differences - quick look only; dic_icgn.m does it properly)
[dux_dx, dux_dy] = gradient(displacement_x, grid_step);
[duy_dx, duy_dy] = gradient(displacement_y, grid_step);
epsilon_xx = dux_dx;
epsilon_yy = duy_dy;
epsilon_xy = 0.5 * (dux_dy + duy_dx);

fprintf('Mean displacement: [dx, dy] = [%.2f, %.2f] px\n', ...
    mean(displacement_x(:),'omitnan'), mean(displacement_y(:),'omitnan'));
fprintf('Mean NCC: %.4f\n', mean(ncc_values(:),'omitnan'));

% Displacement and strain fields
valid_mask = ~isnan(ncc_values);
fields = {displacement_x, displacement_y, ncc_values, epsilon_xx, epsilon_yy, epsilon_xy};
titles = {'Displacement dx (X)','Displacement dy (Y)','NCC', ...
          'Strain \epsilon_{xx}','Strain \epsilon_{yy}','Strain \epsilon_{xy}'};
figure('Name','DIC Results','Position',[100,100,1200,800]);
for k = 1:6
    subplot(2,3,k);
    hImg = imagesc(gx - pad, gy - pad, fields{k});
    set(hImg, 'AlphaData', valid_mask);
    colorbar; axis image;
    title(titles{k}); xlabel('X'); ylabel('Y');
end

% Displacement vector field
figure;
imagesc(def_image); colormap(gray); axis image; hold on;
[Xg, Yg] = meshgrid(gx - pad, gy - pad);
vi = find(valid_mask);
quiver(Xg(vi), Yg(vi), displacement_x(vi), displacement_y(vi), 2, 'r', 'LineWidth', 1.2);
title('Displacement vector field'); xlabel('X'); ylabel('Y'); hold off;

% ===== Local functions =====
function g = to_gray(img)
    if ndims(img) == 3 && size(img,3) >= 3   % >= 3: RGBA PNGs work too
        img = double(img);
        g = 0.2989*img(:,:,1) + 0.5870*img(:,:,2) + 0.1140*img(:,:,3);
    else
        g = double(img);
    end
end

function P = sym_pad(A, p)
    [r, c] = size(A);
    ri = [p:-1:1, 1:r, r:-1:r-p+1];
    ci = [p:-1:1, 1:c, c:-1:c-p+1];
    P = A(ri, ci);
end

function g = read_gray(pathname)
% Grayscale double in [0..1] for any bit depth: integer classes scaled by
% their OWN intmax (uint8 -> /255, uint16 -> /65535, ...); floats above 1
% scaled by their max. Same as /255 for 8-bit input.
    I = imread(pathname);
    if isinteger(I), s = double(intmax(class(I))); else, s = 1; end
    g = to_gray(I) / s;
    mx = max(g(:));
    if mx > 1, g = g / mx; end
end

function names = sort_frames(names)
% Chronological sort: by the trailing number before the extension when
% every name has one (img_2 < img_10), lexicographic otherwise. Zero-padded
% series sort identically either way.
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
