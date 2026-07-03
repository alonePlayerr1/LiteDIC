% Validation of the LiteDIC subpixel engine against DICe (Sandia).
% ONE BUTTON: run this script - everything else is automatic.
%
% The 6_DIC image series IS DICe's official `full_field` example
% (full_field_ref.tiff / full_field_def.tiff are byte-identical to
% plate_ref.tiff / plate_def_01.tiff). DICe ships a gold solution file for
% that pair, bundled here in dice_reference/DICe_solution_0.txt, so the
% first frame is validated against Sandia's official reference numbers
% OUT OF THE BOX - no DICe installation required.
%
% What the script does:
%   1) Runs 6_DIC/dic_icgn.m if its saved results are missing or stale.
%   2) Compares frame 01 against the bundled official gold file.
%   3) If you have run DICe yourself on the full series (dice_results/
%      exists), every available frame is compared as well.
%      To produce those: install DICe (github.com/dicengine/dice), then
%          cd 9_Validation_DICe/dice_input && dice -i input.xml -v
%      The provided input files replicate the official example settings.
%
% Reference settings (from the official example, see dice_reference/):
%   subset 27, step 10, ZNSSD, KEYS_FOURTH interpolation, VSG strain 50 px.
% LiteDIC uses subset 41 and spline interpolation, so small differences are
% expected; agreement criterion for displacements:
%   |bias| < 0.02 px and RMS < 0.05 px per frame.
%
% DICe output format (verified against the gold file): '***' comment lines,
% then a CSV header (SUBSET_ID,COORDINATE_X,...). Columns are located by
% NAME. Both codes use image coordinates (x right, y down). Failed points
% (SIGMA = -1) are excluded.

close all; clc; clear;

%% ============== LiteDIC results (auto-generate if needed) ==============
this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
results_mat = fullfile(this_dir, '..', '6_DIC', 'dic_results.mat');

need_run = ~isfile(results_mat);
if ~need_run
    S = load(results_mat);
    need_run = ~isfield(S, 'series') || ~isfield(S.series, 'U');
    if need_run
        fprintf('dic_results.mat is from an older dic_icgn.m - regenerating.\n');
    end
end
if need_run
    fprintf('Running 6_DIC/dic_icgn.m (takes a moment)...\n');
    % NOTE: dic_icgn.m starts with `clear` - all variables defined above
    % are wiped, so the paths are recomputed from mfilename afterwards.
    evalc('run(fullfile(fileparts(mfilename(''fullpath'')), ''..'', ''6_DIC'', ''dic_icgn.m''))');
    close all;
    this_dir = fileparts(mfilename('fullpath'));
    if isempty(this_dir), this_dir = pwd; end
    results_mat = fullfile(this_dir, '..', '6_DIC', 'dic_results.mat');
    S = load(results_mat);
end
gx = S.R.gx; gy = S.R.gy;
fprintf('LiteDIC results: %d frames, %dx%d grid.\n\n', ...
    numel(S.series), numel(gy), numel(gx));

%% ============== Collect DICe files to compare against ==============
% Each entry: DICe solution file + index into S.series. DICe numbers its
% solution files by deformed-frame index starting at 0, so solution N
% corresponds to series(N+1). Mapping is by the NUMBER in the file name,
% not by sort order (DICe does not zero-pad: _2 sorts after _10).
cmp = struct('file', {}, 'k', {}, 'label', {});

gold = fullfile(this_dir, 'dice_reference', 'DICe_solution_0.txt');
if isfile(gold)
    cmp(end+1) = struct('file', gold, 'k', 1, 'label', 'official gold');
else
    warning('Bundled gold file missing: %s', gold);
end

d = dir(fullfile(this_dir, 'dice_results', 'DICe_solution_*.txt'));
for i = 1:numel(d)
    tok = regexp(d(i).name, 'DICe_solution_(\d+)\.txt$', 'tokens', 'once');
    if isempty(tok), continue; end
    k = str2double(tok{1}) + 1;
    if k < 1 || k > numel(S.series)
        warning('%s has no matching LiteDIC frame (series has %d) - skipped.', ...
            d(i).name, numel(S.series));
        continue;
    end
    cmp(end+1) = struct('file', fullfile(d(i).folder, d(i).name), ...
        'k', k, 'label', 'user-run DICe'); %#ok<SAGROW>
end
if isempty(cmp), error('Nothing to compare against.'); end
if numel(cmp) == 1
    fprintf(['Only the bundled gold frame is available. For the full series,\n' ...
             'run DICe yourself:  cd dice_input && dice -i input.xml -v\n\n']);
end

%% ===================== Per-frame comparison =====================
BIAS_TOL = 0.02;   % px
RMS_TOL  = 0.05;   % px

stats = struct('k', {}, 'frame', {}, 'label', {}, 'npts', {}, ...
               'bias_u', {}, 'bias_v', {}, 'rms_u', {}, 'rms_v', {}, 'pass', {});
fprintf('%-16s %-15s %6s %9s %9s %9s %9s   %s\n', ...
    'LiteDIC frame', 'DICe source', 'pts', 'bias u', 'bias v', 'RMS u', 'RMS v', 'verdict');
for i = 1:numel(cmp)
    k = cmp(i).k;
    [xd, yd, ud, vd] = read_dice_solution(cmp(i).file);

    Fu = griddedInterpolant({gy, gx}, S.series(k).U, 'linear', 'none');
    Fv = griddedInterpolant({gy, gx}, S.series(k).V, 'linear', 'none');
    ul = Fu(yd, xd);   % LiteDIC displacements at the DICe subset centers
    vl = Fv(yd, xd);

    m = isfinite(ul) & isfinite(vl) & isfinite(ud) & isfinite(vd);
    if nnz(m) < 10
        warning('%s: only %d common valid points - skipped.', cmp(i).file, nnz(m));
        continue;
    end
    du = ul(m) - ud(m);
    dv = vl(m) - vd(m);

    % sign-convention guard: the v fields must be positively correlated
    if nnz(m) > 50 && std(vd(m)) > 1e-3 && corr_coeff(vl(m), vd(m)) < 0
        warning(['%s: v fields are ANTI-correlated - y-axis sign ' ...
                 'convention mismatch?'], cmp(i).file);
    end

    st = struct('k', k, 'frame', S.series(k).name, 'label', cmp(i).label, ...
        'npts', nnz(m), 'bias_u', mean(du), 'bias_v', mean(dv), ...
        'rms_u', sqrt(mean(du.^2)), 'rms_v', sqrt(mean(dv.^2)), 'pass', false);
    st.pass = max(abs([st.bias_u, st.bias_v])) < BIAS_TOL && ...
              max([st.rms_u, st.rms_v]) < RMS_TOL;
    stats(end+1) = st; %#ok<SAGROW>
    if st.pass, verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('%-16s %-15s %6d %+9.4f %+9.4f %9.4f %9.4f   %s\n', ...
        st.frame, st.label, st.npts, st.bias_u, st.bias_v, st.rms_u, st.rms_v, verdict);
end
if isempty(stats), error('No frames could be compared.'); end

fprintf(['\nCriterion: |bias| < %.2f px and RMS < %.2f px.\n' ...
         'Overall: %d/%d frames pass. Max |bias| %.4f px, max RMS %.4f px.\n'], ...
    BIAS_TOL, RMS_TOL, nnz([stats.pass]), numel(stats), ...
    max(abs([[stats.bias_u], [stats.bias_v]])), max([[stats.rms_u], [stats.rms_v]]));

%% ===================== Plots =====================
kk = 1:numel(stats);

figure('Name', 'LiteDIC vs DICe: agreement');
subplot(2,1,1);
plot(kk, [stats.bias_u], 'o-', kk, [stats.bias_v], 's-', 'LineWidth', 1.4);
yline( BIAS_TOL, ':r'); yline(-BIAS_TOL, ':r'); yline(0, ':');
grid on; legend('bias u', 'bias v', 'Location', 'best');
ylabel('bias, px'); title('LiteDIC - DICe: bias per compared frame');
subplot(2,1,2);
plot(kk, [stats.rms_u], 'o-', kk, [stats.rms_v], 's-', 'LineWidth', 1.4);
yline(RMS_TOL, ':r'); grid on;
legend('RMS u', 'RMS v', 'Location', 'best');
ylabel('RMS, px'); xlabel('Compared frame');
set(findobj(gcf, 'Type', 'axes'), 'XTick', kk, ...
    'XTickLabel', strrep({stats.frame}, '_', '\_'));

% scatter for the first compared frame (the gold one when present)
[xd, yd, ud, vd] = read_dice_solution(cmp(1).file); %#ok<ASGLU>
Fv = griddedInterpolant({gy, gx}, S.series(stats(1).k).V, 'linear', 'none');
vl = Fv(yd, xd);
m = isfinite(vl) & isfinite(vd);
figure('Name', 'LiteDIC vs DICe: point-by-point scatter');
plot(vd(m), vl(m), '.', 'MarkerSize', 4); hold on;
lims = [min(vd(m)), max(vd(m))];
plot(lims, lims, 'r-', 'LineWidth', 1.2); grid on; axis equal;
xlabel('DICe v, px'); ylabel('LiteDIC v, px');
title(sprintf('v displacement, %s vs %s (line = perfect agreement)', ...
    stats(1).frame, stats(1).label), 'Interpreter', 'none');

save(fullfile(this_dir, 'validation_results.mat'), 'stats');
fprintf('Saved per-frame statistics to validation_results.mat\n');

%% ===================== Local functions =====================

function [x, y, u, v] = read_dice_solution(fname)
% Read one DICe solution file. Format verified against the official gold
% file: '***' comment lines, then a CSV header row; columns are found by
% NAME so order or extra fields do not matter. Points flagged invalid by
% DICe (SIGMA = -1) are returned as NaN.
    fid = fopen(fname, 'r');
    if fid < 0, error('Cannot open %s', fname); end
    n_skip = 0; header = '';
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        t = strtrim(line);
        if isempty(t) || startsWith(t, '*') || startsWith(t, '%') || startsWith(t, '#')
            n_skip = n_skip + 1;
            continue;
        end
        header = t;       % first non-comment line = column names
        break;
    end
    fclose(fid);
    if isempty(header), error('No header row found in %s', fname); end

    cols = strtrim(strsplit(header, ','));
    ix = find_col(cols, 'COORDINATE_X', fname);
    iy = find_col(cols, 'COORDINATE_Y', fname);
    iu = find_col(cols, 'DISPLACEMENT_X', fname);
    iv = find_col(cols, 'DISPLACEMENT_Y', fname);
    is = find(strcmpi(cols, 'SIGMA'), 1);   % optional

    M = readmatrix(fname, 'NumHeaderLines', n_skip + 1, ...
        'Delimiter', ',', 'FileType', 'text');
    if size(M, 2) < numel(cols)
        error('Data/header column mismatch in %s', fname);
    end
    x = M(:, ix); y = M(:, iy); u = M(:, iu); v = M(:, iv);
    if ~isempty(is)
        bad = M(:, is) < 0;                 % DICe marks failed points sigma=-1
        u(bad) = NaN; v(bad) = NaN;
    end
end

function i = find_col(cols, name, fname)
    i = find(strcmpi(cols, name), 1);
    if isempty(i)
        error('Column %s not found in %s (header: %s)', ...
            name, fname, strjoin(cols, ','));
    end
end

function c = corr_coeff(a, b)
% Pearson correlation without toolboxes.
    a = a(:) - mean(a); b = b(:) - mean(b);
    c = (a.' * b) / max(norm(a) * norm(b), eps);
end
