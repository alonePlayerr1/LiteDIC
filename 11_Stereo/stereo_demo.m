% Stereo (3D) DIC demonstration on the two-camera application pair
% 9_Validation_DICe/stereo_sample/ (GT4-0110 -> GT4-0129). ONE BUTTON.
%
% Unlike stereo_d_spec (see validate_dice_stereo.m) this pair has NO
% reference solution, so no accuracy is claimed here: the run is judged
% by internal consistency only - epipolar residuals of both matched
% pairs under the calibrated geometry, match quality, and the physical
% plausibility of the reconstructed surface and displacement field.
% GT4-0110 is used as the reference state; the frames come from the same
% experiment as stereo_d_spec, ~105 load steps later.
%
% Outputs: printed summary + stereo_demo_results.mat +
% stereo_demo_surface.png / stereo_demo_d3d.png in 9_Validation_DICe/.

close all; clc; clear;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root = fileparts(this_dir);
addpath(this_dir, fullfile(root, '6_DIC'));
sdir = fullfile(root, '9_Validation_DICe', 'stereo_sample');
odir = fullfile(root, '9_Validation_DICe');

%% ===================== geometry + parameters =====================
G = stereo_geometry();
fprintf('Geometry: baseline %.2f mm | pitch %.3f mm | toe-in %.2f deg\n\n', ...
    G.baseline_mm, G.pitch_mm, G.E.ang_deg);

P = struct();
P.subset_size     = 25;
P.grid_step       = 20;
P.search_range    = 40;      % same regime as the validated stereo_d_spec run
P.search_temporal = 25;      % 19 load steps between the frames - allow more
P.black_thr       = 0.15;
P.zncc_min        = 0.75;
P.zncc_min_cross  = 0.50;
P.max_iter        = 30;
P.conv_tol        = 1e-3;

%% ===================== run =====================
f0 = {fullfile(sdir, 'GT4-0110_0.tif'), fullfile(sdir, 'GT4-0129_0.tif')};
f1 = {fullfile(sdir, 'GT4-0110_1.tif'), fullfile(sdir, 'GT4-0129_1.tif')};
t0 = tic;
S = stereo_dic(f0, f1, G, P);
fprintf('stereo_dic: 2 states in %.1f s\n\n', toc(t0));

%% ===================== internal-consistency summary =====================
for k = 1:2
    D = S.frames(k).D;
    q = D.Q(isfinite(D.Q));
    fprintf(['state %d: shift [%4d %4d] | valid %5d | cross ZNCC median ' ...
        '%.3f | epi med %.3f px (RMS %.3f)\n'], k-1, D.shift(1), D.shift(2), ...
        S.frames(k).nvalid, median(q), S.frames(k).epi_med, S.frames(k).epi_rms);
end

U0 = S.frames(2).U0; V0 = S.frames(2).V0;
fprintf('temporal cam0 (0110 -> 0129): median u %+.3f px, v %+.3f px, valid %d\n', ...
    median(U0(isfinite(U0))), median(V0(isfinite(V0))), nnz(isfinite(U0)));

d3 = S.frames(2).d3D_mm;
mag = sqrt(sum(d3.^2, 3));
dz = d3(:, :, 3);
fprintf(['3D displacement: median |d3D| %.4f mm | mean dZ %+.4f mm | ' ...
    'valid %d\n'], median(mag(isfinite(mag))), mean(dz(isfinite(dz))), ...
    nnz(isfinite(mag)));

Z = S.X_ref_mm(:, :, 3);
fprintf('reference surface: Z range %.2f .. %.2f mm (%.1f%% of nodes valid)\n\n', ...
    min(Z(:), [], 'omitnan'), max(Z(:), [], 'omitnan'), ...
    100 * nnz(isfinite(Z)) / numel(Z));

%% ===================== figures =====================
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 720 540]);
surf(S.X_ref_mm(:,:,1), S.X_ref_mm(:,:,2), S.X_ref_mm(:,:,3), ...
    S.X_ref_mm(:,:,3), 'EdgeColor', 'none');
axis equal tight; view(-35, 25); colorbar;
xlabel('X, mm'); ylabel('Y, mm'); zlabel('Z, mm');
title('Reference surface (GT4-0110), camera-0 frame');
print(fig, fullfile(odir, 'stereo_demo_surface.png'), '-dpng', '-r150');

fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 720 540]);
imagesc(S.gx, S.gy, mag, 'AlphaData', isfinite(mag));
axis image; set(gca, 'YDir', 'reverse'); colorbar;
xlabel('x, px (cam0)'); ylabel('y, px (cam0)');
title('|d3D|, mm (GT4-0110 \rightarrow GT4-0129)');
print(fig2, fullfile(odir, 'stereo_demo_d3d.png'), '-dpng', '-r150');

save(fullfile(odir, 'stereo_demo_results.mat'), 'S', 'P', 'G');
fprintf(['Saved stereo_demo_results.mat, stereo_demo_surface.png, ' ...
    'stereo_demo_d3d.png.\nNo reference solution exists for this pair: ' ...
    'the numbers above are internal-consistency\nindicators, not accuracy ' ...
    'claims (see validate_dice_stereo.m for the gold-checked run).\n']);
