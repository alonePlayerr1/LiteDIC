function R = validate_external(src, opts)
% VALIDATE_EXTERNAL  Audit a third-party DIC displacement field against
% the DICe gold standard bundled with this repository.
%
%   R = VALIDATE_EXTERNAL(src) takes the results of ANY external DIC code
%   (src = path to a CSV/TSV file, or the struct from IMPORT_EXTERNAL_DIC)
%   computed on the bundled DICe example (9_Validation_DICe, reference
%   frame -> deformed frame 01) and compares them point-by-point with the
%   official DICe gold solution:
%     * the external field is interpolated linearly (NO extrapolation)
%       onto the 972 gold subset centres;
%     * bias and RMS of u and v are computed exactly as in
%       VALIDATE_DICE (the same PASS gates: |bias| < 0.02 px,
%       RMS < 0.05 px);
%     * an error map figure and a plain-text summary are written next to
%       the input (or to opts.outdir).
%
%   Any DIC code (classical or neural) can be checked this way: run it on
%   the bundled images and compare against the DICe gold (Sandia).
%
%   opts (all optional):
%     .outdir     folder for report files (default: folder of src / pwd)
%     .cols       [ix iy iu iv] for the importer
%     .bias_max   PASS gate for |bias|, px (0.02)
%     .rms_max    PASS gate for RMS, px (0.05)
%     .gold       path to a DICe solution file (default: bundled gold 0)
%
%   Output struct R: n_common, bias_u/v, rms_u/v, med_abs_u/v, pass,
%   report (text), files (written paths).
%
%   See also IMPORT_EXTERNAL_DIC, VALIDATE_DICE.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    if nargin < 2, opts = struct(); end
    bias_max = getf(opts, 'bias_max', 0.02);
    rms_max  = getf(opts, 'rms_max', 0.05);
    goldfile = getf(opts, 'gold', fullfile(root, '9_Validation_DICe', ...
        'dice_reference', 'DICe_solution_0.txt'));

    % ---- load both fields ----
    if ischar(src) || isstring(src)
        E = import_external_dic(char(src), opts);
        srcname = char(src);
    elseif isstruct(src)
        E = src; srcname = getf(getf(E, 'info', struct()), 'file', 'in-memory field');
    else
        error('validate_external:src', 'src must be a file path or an imported struct.');
    end
    G = import_external_dic(goldfile);          % gold parses with the same reader
    if isempty(getf(opts, 'outdir', ''))
        outdir = fileparts(srcname);
        if isempty(outdir), outdir = pwd; end
    else
        outdir = opts.outdir; if ~isfolder(outdir), mkdir(outdir); end
    end

    % ---- interpolate the external field onto the gold points ----
    m = isfinite(E.u) & isfinite(E.v);
    if nnz(m) < 8
        error('validate_external:few', 'Only %d valid external points.', nnz(m));
    end
    Fu = scatteredInterpolant(E.x(m), E.y(m), E.u(m), 'linear', 'none');
    Fv = scatteredInterpolant(E.x(m), E.y(m), E.v(m), 'linear', 'none');
    ug = Fu(G.x, G.y); vg = Fv(G.x, G.y);

    ok = isfinite(ug) & isfinite(vg) & isfinite(G.u) & isfinite(G.v);
    du = ug(ok) - G.u(ok); dv = vg(ok) - G.v(ok);
    R.n_common  = nnz(ok);
    R.bias_u    = mean(du);          R.bias_v    = mean(dv);
    R.rms_u     = sqrt(mean(du.^2)); R.rms_v     = sqrt(mean(dv.^2));
    R.med_abs_u = median(abs(du));   R.med_abs_v = median(abs(dv));
    R.pass = abs(R.bias_u) < bias_max && abs(R.bias_v) < bias_max && ...
             R.rms_u < rms_max && R.rms_v < rms_max && R.n_common >= 100;

    % ---- error map figure ----
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 1000 420]);
    subplot(1, 2, 1);
    scatter(G.x(ok), G.y(ok), 14, du, 'filled');
    axis image; set(gca, 'YDir', 'reverse'); colorbar;
    title('u error vs DICe gold, px'); xlabel('x, px'); ylabel('y, px');
    subplot(1, 2, 2);
    scatter(G.x(ok), G.y(ok), 14, dv, 'filled');
    axis image; set(gca, 'YDir', 'reverse'); colorbar;
    title('v error vs DICe gold, px'); xlabel('x, px'); ylabel('y, px');
    f_png = fullfile(outdir, 'external_error_map.png');
    print(fig, f_png, '-dpng', '-r150'); close(fig);

    % ---- text report ----
    if R.pass, verdict = 'PASS'; else, verdict = 'FAIL'; end
    L = {sprintf('EXTERNAL DIC VALIDATION vs DICe gold - %s', char(datetime('now')))
         sprintf('input : %s', srcname)
         sprintf('gold  : %s', goldfile)
         sprintf('common points: %d', R.n_common)
         sprintf('bias u: %+.4f px | bias v: %+.4f px  (gate |bias| < %.3g)', ...
                 R.bias_u, R.bias_v, bias_max)
         sprintf('RMS  u: %.4f px  | RMS  v: %.4f px  (gate RMS < %.3g)', ...
                 R.rms_u, R.rms_v, rms_max)
         sprintf('median |err|: u %.4f px, v %.4f px', R.med_abs_u, R.med_abs_v)
         sprintf('VERDICT: %s', verdict)};
    R.report = strjoin(L, newline);
    f_txt = fullfile(outdir, 'external_validation.txt');
    fid = fopen(f_txt, 'w'); fwrite(fid, R.report); fclose(fid);
    R.files = {f_png, f_txt};
    fprintf('%s\n', R.report);
    fprintf('Report: %s\n        %s\n', f_txt, f_png);
end

function v = getf(s, name, dflt)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name);
    else, v = dflt; end
end
