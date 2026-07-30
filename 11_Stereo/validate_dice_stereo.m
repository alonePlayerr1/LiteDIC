% Validation of the LiteDIC stereo pipeline against DICe's own stereo
% reference results (gold) on the GT4 series. ONE BUTTON: run this script.
%
% Data: 9_Validation_DICe/stereo_d_spec/ (from the DICe repository,
% tests/regression/stereo_d_spec): synchronized pairs GT4-0000..0005_{0,1},
% DICe settings (subset 25, step 100, KEYS_FOURTH), and DICe's gold
% solutions gold/DICe_solution_0..5.txt. Gold columns (located by NAME,
% '***' comments skipped, failed points MATCH = -1 excluded):
%   COORDINATE_X/Y     cam0 subset centre (px, image coords, x right y down)
%   CROSS_CORR_Q/R     cam0 -> cam1 cross-camera match offset (px)
%   DISPLACEMENT_X/Y   temporal displacement of the cam0 subset (px, vs frame 0)
%
% What is compared, per frame:
%   1) disparity  : LiteDIC cross-camera field, interpolated at the gold
%                   centres, vs (CROSS_CORR_Q, CROSS_CORR_R) - bias/RMS.
%   2) temporal   : LiteDIC cam0 displacement vs DISPLACEMENT_X/Y - bias/RMS.
%   3) epipolar   : RMS distance to the epipolar lines of OUR calibrated
%                   geometry, for gold pairs and for LiteDIC pairs. This
%                   simultaneously checks the geometry and the Q/R sign
%                   convention: a wrong convention explodes this number.
%   4) 3D surface : both correspondence sets triangulated with the SAME
%                   geometry; RMS difference of the surfaces (mm). Using
%                   one geometry for both isolates matching differences
%                   from calibration.
%
% Numbers are PRINTED (and saved to
% 9_Validation_DICe/stereo_validation_results.mat); RUN_TESTS intervals
% get pinned to them afterwards, per the repo convention. Nothing here
% modifies the 2D engine.

close all; clc; clear;

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
root = fileparts(this_dir);
addpath(this_dir, fullfile(root, '6_DIC'));

spec = fullfile(root, '9_Validation_DICe', 'stereo_d_spec');
nframes = 6;                                   % GT4-0000 .. GT4-0005

%% ===================== geometry + DIC parameters =====================
G = stereo_geometry();
fprintf(['Geometry: baseline %.2f mm | pitch %.3f mm | toe-in %.2f deg | ' ...
    'common views %d\n\n'], G.baseline_mm, G.pitch_mm, G.E.ang_deg, G.E.nviews);

P = struct();
P.subset_size     = 25;      % = DICe gold subset size
P.grid_step       = 20;      % dense enough to interpolate at step-100 centres
P.search_range    = 40;      % residual after the global stereo pre-shift.
                             % Residuals up to ~35 px
                             % (global FFT shift [-22 -42] vs disparities
                             % -46..-68), so 25 left the seed search at the
                             % basin edge; 40 gives margin, and the gold
                             % centres (y >= 124, x in 124..1824) stay
                             % inside the engine margin (12+40+2 = 54).
P.search_temporal = 10;      % temporal motion in this series is small
P.black_thr       = 0.15;    % masks the dark surround (max-in-window test)
P.zncc_min        = 0.75;    % temporal quality gate
P.zncc_min_cross  = 0.50;    % cross-camera: perspective lowers ZNCC
P.max_iter        = 30;
P.conv_tol        = 1e-3;

%% ===================== run the stereo pipeline =====================
f0 = cell(1, nframes); f1 = cell(1, nframes);
for k = 1:nframes
    f0{k} = fullfile(spec, 'images', sprintf('GT4-%04d_0.tif', k-1));
    f1{k} = fullfile(spec, 'images', sprintf('GT4-%04d_1.tif', k-1));
end
t0 = tic;
S = stereo_dic(f0, f1, G, P);
fprintf('stereo_dic: %d frames in %.1f s\n', nframes, toc(t0));
for k = 1:nframes
    fprintf(['  frame %d: shift [%4d %4d] | valid %5d | epi RMS %.3f ' ...
        'med %.3f px | peakq %.2f\n'], k-1, S.frames(k).D.shift(1), ...
        S.frames(k).D.shift(2), S.frames(k).nvalid, S.frames(k).epi_rms, ...
        S.frames(k).epi_med, S.frames(k).D.peakq);
end
fprintf('\n');

%% ===================== per-frame comparison vs gold =====================
res = struct('frame', {}, 'npts', {}, ...
    'bias_qx', {}, 'bias_qy', {}, 'rms_qx', {}, 'rms_qy', {}, 'med_q', {}, ...
    'bias_u', {}, 'bias_v', {}, 'rms_u', {}, 'rms_v', {}, ...
    'epi_gold', {}, 'epi_ours', {}, 'rms3d_mm', {}, 'dz_mm', {});

% NOTE on frames 1..5: gold's CROSS_CORR is FRAME-CONSTANT (DICe matches
% the cameras once, on the reference pair - verified in the diagnostics
% below), while LiteDIC re-matches every synchronized pair. Any real
% motion of the specimen therefore appears as a growing LiteDIC-vs-gold
% disparity "bias" that is NOT an error; the dZ column (our mean
% out-of-plane displacement) is printed to show the two track together.
% The apples-to-apples disparity comparison is frame 0.
fprintf('%-6s %5s | %-31s | %-25s | %-13s | %-9s | %s\n', 'frame', 'pts', ...
    'disparity bias/RMS x,y / med (px)', 'temporal bias/RMS x,y (px)', ...
    'epi gold/ours', '3D RMS mm', 'dZ mm');
for k = 1:nframes
    gfile = fullfile(spec, 'gold', sprintf('DICe_solution_%d.txt', k-1));
    T = read_dice_gold(gfile);
    v = T.MATCH >= 0;                       % gold's own failed points out
    xg = T.COORDINATE_X(v); yg = T.COORDINATE_Y(v);
    qg = T.CROSS_CORR_Q(v); rg = T.CROSS_CORR_R(v);
    ug = T.DISPLACEMENT_X(v); vg = T.DISPLACEMENT_Y(v);

    % LiteDIC fields at the gold centres
    Dk = S.frames(k).D;
    Fu = griddedInterpolant({Dk.gy, Dk.gx}, Dk.U, 'linear', 'none');
    Fv = griddedInterpolant({Dk.gy, Dk.gx}, Dk.V, 'linear', 'none');
    ql = Fu(yg, xg); rl = Fv(yg, xg);

    Tu = griddedInterpolant({S.gy, S.gx}, S.frames(k).U0, 'linear', 'none');
    Tv = griddedInterpolant({S.gy, S.gx}, S.frames(k).V0, 'linear', 'none');
    ul = Tu(yg, xg); vl = Tv(yg, xg);

    m = isfinite(ql) & isfinite(rl);
    md = m & isfinite(ul) & isfinite(vl);

    % 1) disparity
    dqx = ql(m) - qg(m); dqy = rl(m) - rg(m);
    % 2) temporal
    dux = ul(md) - ug(md); dvy = vl(md) - vg(md);
    % 3) epipolar, same geometry for both
    epi_g = epi_rms([xg(m), yg(m)], [xg(m)+qg(m), yg(m)+rg(m)], G.F);
    epi_l = epi_rms([xg(m), yg(m)], [xg(m)+ql(m), yg(m)+rl(m)], G.F);
    % 4) 3D surfaces, same geometry for both
    Xl = triangulate_dlt([xg(m), yg(m)], [xg(m)+ql(m), yg(m)+rl(m)], ...
        G.K0, G.d0, G.K1, G.d1, G.R, G.T) * G.pitch_mm;
    Xd = triangulate_dlt([xg(m), yg(m)], [xg(m)+qg(m), yg(m)+rg(m)], ...
        G.K0, G.d0, G.K1, G.d1, G.R, G.T) * G.pitch_mm;
    rms3d = sqrt(mean(sum((Xl - Xd).^2, 2)));

    med_q = median(hypot(dqx, dqy));
    dz = S.frames(k).d3D_mm(:, :, 3);
    dz = mean(dz(isfinite(dz)));

    res(k).frame = k-1;              res(k).npts = nnz(m);
    res(k).bias_qx = mean(dqx);      res(k).bias_qy = mean(dqy);
    res(k).rms_qx = rms(dqx);        res(k).rms_qy = rms(dqy);
    res(k).med_q = med_q;
    res(k).bias_u = mean(dux);       res(k).bias_v = mean(dvy);
    res(k).rms_u = rms(dux);         res(k).rms_v = rms(dvy);
    res(k).epi_gold = epi_g;         res(k).epi_ours = epi_l;
    res(k).rms3d_mm = rms3d;         res(k).dz_mm = dz;

    fprintf(['%-6d %5d | %+7.3f %+7.3f %6.3f %6.3f %6.3f | ' ...
        '%+7.4f %+7.4f %6.4f %6.4f | %5.3f %5.3f | %8.5f | %+8.5f\n'], ...
        k-1, nnz(m), mean(dqx), mean(dqy), rms(dqx), rms(dqy), med_q, ...
        mean(dux), mean(dvy), rms(dux), rms(dvy), epi_g, epi_l, rms3d, dz);
end

%% ============== frame-0 diagnostics (RMS 0.356 px there ===============
% vs 0.09-0.20 elsewhere; find out whether that is outliers or field-wide,
% and whether gold's CROSS_CORR is frame-constant).
fprintf('\n--- frame-0 diagnostics ---\n');
T0 = read_dice_gold(fullfile(spec, 'gold', 'DICe_solution_0.txt'));
T1 = read_dice_gold(fullfile(spec, 'gold', 'DICe_solution_1.txt'));
[c01, i0, i1] = intersect(T0.SUBSET_ID, T1.SUBSET_ID);
dQ = T0.CROSS_CORR_Q(i0) - T1.CROSS_CORR_Q(i1);
dR = T0.CROSS_CORR_R(i0) - T1.CROSS_CORR_R(i1);
fprintf('gold cross frame0-vs-frame1: max |dQ| %.4f px, max |dR| %.4f px %s\n', ...
    max(abs(dQ)), max(abs(dR)), tern(max(abs([dQ;dR])) < 1e-6, ...
    '(FRAME-CONSTANT - DICe cross-matches once)', '(varies per frame)'));

v0 = T0.MATCH >= 0;
xg = T0.COORDINATE_X(v0); yg = T0.COORDINATE_Y(v0);
qg = T0.CROSS_CORR_Q(v0); rg = T0.CROSS_CORR_R(v0);
qrl = cell(1, 2);                       % our frames 0 and 1 at gold centres
for k = 1:2
    Dk = S.frames(k).D;
    Fu = griddedInterpolant({Dk.gy, Dk.gx}, Dk.U, 'linear', 'none');
    Fv = griddedInterpolant({Dk.gy, Dk.gx}, Dk.V, 'linear', 'none');
    qrl{k} = [Fu(yg, xg), Fv(yg, xg)];
end
m01 = all(isfinite([qrl{1}, qrl{2}]), 2);
d01 = qrl{1}(m01,:) - qrl{2}(m01,:);
fprintf('OUR frame0-vs-frame1 disparity at gold centres: RMS x %.3f  y %.3f px\n', ...
    rms(d01(:,1)), rms(d01(:,2)));

ex = qrl{1}(:,1) - qg;  ey = qrl{1}(:,2) - rg;
mm0 = isfinite(ex) & isfinite(ey);
r2 = hypot(ex(mm0), ey(mm0));
madr = 1.4826 * median(abs(r2 - median(r2)));
fprintf(['frame-0 residual vs gold: median %.3f px | robust sigma %.3f | ' ...
    'max %.3f px | pts > 3*sigma: %d of %d\n'], ...
    median(r2), madr, max(r2), nnz(r2 > median(r2) + 3*madr), numel(r2));
[~, order] = sort(r2, 'descend');
xv = xg(mm0); yv = yg(mm0);
fprintf('worst 5 residuals (x, y, |d|):\n');
for i = 1:min(5, numel(order))
    j = order(i);
    fprintf('  (%6.0f, %5.0f)  %.3f px\n', xv(j), yv(j), r2(j));
end

save(fullfile(root, '9_Validation_DICe', 'stereo_validation_results.mat'), ...
    'res', 'P', 'G', 'S');
fprintf(['\nSaved stereo_validation_results.mat. RUN_TESTS intervals are\n' ...
    'set from these values.\n' ...
    'Interpretation aids: epi gold vs ours should be SIMILAR (both limited\n' ...
    'by the calibration floor ~0.5 px); a huge gold epipolar RMS would mean\n' ...
    'the CROSS_CORR_Q/R sign convention is wrong.\n']);

%% ===================== local functions =====================
function r = epi_rms(uv0, uv1, F)
    n = size(uv0, 1);
    l = F * [uv0, ones(n,1)].';
    num = sum([uv1, ones(n,1)].' .* l, 1);
    den = hypot(l(1,:), l(2,:));
    r = sqrt(mean((num ./ den).^2));
end

function y = rms(x)
    y = sqrt(mean(x.^2));
end

function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end

