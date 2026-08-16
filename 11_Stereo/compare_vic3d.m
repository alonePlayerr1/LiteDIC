function out = compare_vic3d(dataDir)
% COMPARE_VIC3D  Compare LiteDIC stereo deflection profiles with a Vic-3D
% line export and recover the metric scale of the unit-pitch calibration.
%
%   compare_vic3d('C:\...\DIC 3D dataset')
%   out = compare_vic3d('C:\...\DIC 3D dataset')
%
% With an output argument the same diagnostic is returned as a struct so the
% verdict numbers can be checked programmatically (e.g. by RUN_TESTS):
%   out.s_z         depth agreement (1.0 = perfect), = abs(best.sz)
%   out.s_xy        lateral scale / pitch, mm (= 1 in a metric run)
%   out.b           line-origin offset, mm
%   out.orientation +1 same direction, -1 flipped
%   out.metric      true if the run was already metric (VicCal or Pitch~=1)
%   out.rms         per-state rms of the mm overlay, mm (row vector)
%   out.rms_max     max over states, mm
%   out.drift       per-state whole-rig drift, mm (row vector)
%   out.states      state stems (cell)
%   out.s_z_speckle    depth agreement fitted on SPECKLED rows only
%   out.rms_speckle    per-state overlay rms on speckled rows, mm
%   out.rms_max_speckle max over states, mm
%   out.n_speckle      number of speckled profile rows used
%   out.cover          fraction of the Vic line our measured window spans
%   out.curv_vic       curvature (detrended std) on the window BOTH cover, mm
%   out.curv_ours      the same for our profile - comparable, unlike the raw
%                      table, which detrends each side over its own span
%   out.agreed         |s_z_speckle - 1| <= 0.05 and rms_max_speckle <= 5 um
%   out.pass           the binding verdict: agreed AND cover >= 0.80
%
% Coverage is in the verdict because everything else is computed on the
% overlap of the two lines: half a specimen compared with the same half of
% the reference looks the same as a full one (DIC 3D dataset: s_z 0.993,
% 0.6 um worst residual with half measured).
%
% Speckle-only verdict: Vic-3D (like any DIC) only measures where there is
% speckle, over bare regions its line export is interpolated. Comparing
% against that makes no sense (weld dataset: apparent depth error 15%,
% painted zones agree to ~2 um at s_z 0.98). So the binding verdict uses
% only rows with speckle texture in the reference image; the all-row fit
% is printed for context.
% Report and figures are the same with or without the output argument.
%
% Inputs (produced earlier):
%   <dataDir>\litedic_out\stereo_results.mat   from RUN_STEREO_DATASET
%   <dataDir>\21r.csv                          Vic-3D export: per test frame
%       a "L0" line profile, columns Y [mm]; W [mm] (out-of-plane)
%
% Folder layout: the reference companion files (21r.csv, project.xml) are
% looked up in dataDir first and then in its PARENT folder, so both the
% original flat layout and the split layout used with RUN_STEREO_DATASET's
% 'CalDir' option (<set>\Dataset\ for images + <set>\ for the exports)
% work unchanged.
%
% Method. Everything LiteDIC reports in pitch units scales by ONE unknown
% factor s (the target dot pitch, mm) - the along-beam coordinate and the
% deflection alike:   s * W_ours,k(y) == +/- W_vic,k( s*y + b ).
% The fit is GLOBAL: a single (s, b, orientation) must explain ALL load
% states simultaneously - the growth of the deflection amplitude with
% load pins the scale far harder than any single profile. Before fitting,
% grid rows whose deflection does not follow the common load evolution
% (static fixtures caught by the ROI, clutter) are discarded: per row, the
% correlation of W(row, state) with the median load curve must be >= 0.8.
% Prints the global fit and its per-state residuals (mm), saves an
% absolute-mm overlay to litedic_out\vic3d_compare.png. Diagnostic only.

    outDir = fullfile(dataDir, 'litedic_out');
    L = load(fullfile(outDir, 'stereo_results.mat'));
    B = parse_vic_csv(aux_file(dataDir, '21r.csv'));
    fprintf('Vic-3D export: %d frame blocks\n', numel(B));

    % ---- our data translated into Vic's convention ----
    % Vic-3D (autoplanefit) reports W and Y in the best-fit-plane
    % coordinate system of the reference surface; our raw dZ and camera-Y
    % differ from that by the ~24 deg surface tilt. Translate: fit a plane
    % to our own X_ref, take W = d3D . n (normal component) and
    % Y = in-plane coordinate along the beam (principal direction).
    X = L.S.X_ref_mm;
    P3 = reshape(X, [], 3);
    v = all(isfinite(P3), 2);
    c0 = mean(P3(v, :), 1);
    [~, ~, Vv] = svd(P3(v, :) - c0, 0);
    e1 = Vv(:, 1);                 % along the beam (largest extent)
    n  = Vv(:, 3);                 % surface normal
    fprintf('Reference plane: tilt to camera Z = %.1f deg\n', ...
        acosd(abs(n(3))));
    nY = size(X, 1);
    yo = nan(nY, 1);
    for r = 1:nY
        xr = squeeze(median(X(r, :, :), 2, 'omitnan'));
        yo(r) = (xr(:) - c0(:)).' * e1;
    end
    Wo = []; Vic = {}; names = {};
    for k = 2:numel(L.S.frames)
        stem = L.test_stems{L.sel(k)};
        bi = find(strcmp({B.stem}, [stem '_0']), 1);
        if isempty(bi), continue; end
        D3 = L.S.frames(k).d3D_mm;
        dn = D3(:,:,1) * n(1) + D3(:,:,2) * n(2) + D3(:,:,3) * n(3);
        Wo(:, end+1) = median(dn, 2, 'omitnan'); %#ok<AGROW>
        [yv, iu] = unique(B(bi).Y);
        Vic{end+1} = struct('Y', yv, 'W', B(bi).W(iu)); %#ok<AGROW>
        names{end+1} = stem;                     %#ok<AGROW>
    end
    ns = numel(names);
    if ns < 2, error('compare_vic3d: need >= 2 matched states.'); end
    fprintf('Matched states: %s\n', strjoin(names, ', '));

    % ---- row screen: keep rows that follow the common load evolution ----
    fin = all(isfinite(Wo), 2) & isfinite(yo);
    amp = median(abs(Wo(fin, :)), 1).';                 % load curve
    rc = nan(size(yo));
    for r = find(fin).'
        rc(r) = corr_safe(Wo(r, :).', amp);
    end
    keep = fin & abs(rc) >= 0.8;   % |corr|: our dZ sign convention may be
                                   % opposite to the |W| load curve
    % restrict to the SPECIMEN AOI from the Vic project (the frame also
    % contains the testing machine, which Vic masked with polygons; the
    % first polygon is the specimen gauge region)
    % All polygons, not only the first: a specimen painted in separate zones
    % has one polygon per zone in Vic, reading the first restricts the
    % comparison to that zone (DIC 3D dataset). Rows between polygons stay out:
    % Vic does not measure there either, its export just bridges them.
    xmlf = aux_file(dataDir, 'project.xml');
    if isfile(xmlf)
        ol = regexp(fileread(xmlf), '<outline>([^<]*)</outline>', 'tokens');
        if ~isempty(ol)
            inAoi = false(size(L.S.gy(:)));
            txt = cell(1, numel(ol));
            for q = 1:numel(ol)
                vtx = sscanf(ol{q}{1}, '%f');
                ys = vtx(2:2:end);
                inAoi = inAoi | (L.S.gy(:) >= min(ys) & L.S.gy(:) <= max(ys));
                txt{q} = sprintf('%.0f..%.0f', min(ys), max(ys));
            end
            keep = keep & inAoi;
            fprintf('Specimen AOI from project.xml: %d polygon(s), image rows %s px\n', ...
                numel(ol), strjoin(txt, ' + '));
        end
    end
    fprintf('Rows: %d finite, %d follow the load evolution (kept)\n', ...
        nnz(fin), nnz(keep));
    if nnz(keep) < 10
        error(['compare_vic3d: almost no rows follow the load evolution - ' ...
               'check that the states span different loads.']);
    end
    yo = yo(keep) - median(yo(keep));
    Wo = Wo(keep, :);

    % ---- speckle map: which kept rows lie on painted (measurable) bands ----
    gyk = L.S.gy(:); gyk = gyk(keep);
    spk = speckle_rows(dataDir, L, gyk);
    fprintf('Speckle rows (texture of the reference frame): %d of %d kept\n', ...
        nnz(spk), numel(spk));

    % ---- global coarse-to-fine scan over (s_xy, b, orientation) ----
    % The depth (deflection) scale s_z is solved ANALYTICALLY per
    % candidate by least squares, so the scan is 2D and the fit matches
    % the profile SHAPE, not just an amplitude corner. If the pipeline
    % were perfectly consistent, s_z == s_xy (one physical pitch); their
    % ratio measures the depth/lateral anisotropy of the calibration
    % (a weak-perspective handheld set pins the focal length poorly, and
    % a focal bias stretches depth relative to lateral coordinates).
    % High coverage is required first so the fit spans the WHOLE line;
    % relaxed only if no solution exists.
    % ---- fit-free diagnostic: curvature amplitude of each side ----
    % std of the linearly detrended profile over each line SEPARATELY -
    % no alignment involved, so this is a raw estimate of the depth
    % agreement if both lines span (mostly) the same physical beam.
    fprintf('\nRaw curvature amplitudes (std of detrended profile):\n');
    fprintf('%-10s %10s %10s %8s\n', 'state', 'vic(mm)', 'ours(mm)', 'ratio');
    mfin = isfinite(yo) & all(isfinite(Wo), 2);
    Qo = [ones(nnz(mfin), 1), yo(mfin)];
    for k = 1:ns
        Qv = [ones(numel(Vic{k}.Y), 1), Vic{k}.Y];
        wdv = Vic{k}.W - Qv * (Qv \ Vic{k}.W);
        wdo = Wo(mfin, k) - Qo * (Qo \ Wo(mfin, k));
        fprintf('%-10s %10.4f %10.4f %8.3f\n', names{k}, ...
            std(wdv), std(wdo), std(wdv) / std(wdo));
    end
    span_vic  = max(Vic{1}.Y) - min(Vic{1}.Y);
    span_ours = max(yo(mfin)) - min(yo(mfin));
    fprintf('Line spans: vic %.1f mm | ours %.1f mm\n', span_vic, span_ours);
    % Ratios above detrend each side over its own span, only meaningful when
    % spans agree. DIC 3D dataset: 131 mm vs 61 mm, ratios 23x..105x, but on
    % the shared window Vic's curvature equals ours.
    if abs(span_ours - span_vic) > 0.10 * span_vic
        fprintf(['  the spans differ by %.0f%% - the ratios above compare ' ...
            'different lengths of beam\n  and are NOT a measure of ' ...
            'agreement; see the matched-window table below\n'], ...
            100 * abs(span_ours - span_vic) / span_vic);
    end

    % raw side-by-side view (no fitting at all)
    figr = figure('Visible', 'off', 'Position', [40 40 1100 450]);
    subplot(1, 2, 1); hold on;
    for k = 1:ns, plot(Vic{k}.Y, Vic{k}.W, 'LineWidth', 1.2); end
    grid on; xlabel('Y, mm'); ylabel('W, mm'); title('Vic-3D (raw)');
    subplot(1, 2, 2); hold on;
    for k = 1:ns, plot(yo(mfin), Wo(mfin, k), 'LineWidth', 1.2); end
    grid on; xlabel('Y_{plane}, mm'); ylabel('W = d3D\cdotn, mm');
    title('LiteDIC (raw, plane CS)');
    print(figr, fullfile(outDir, 'vic3d_raw.png'), '-dpng', '-r120');
    close(figr);

    % When the stereo run is already METRIC (Vic calibration or a known
    % pitch), the lateral scale is not a free parameter: s_xy = 1 by
    % construction, and fitting it would let the near-linear profiles
    % slide onto arbitrary subsegments. Only the offset b (line origin)
    % and the orientation remain free.
    metric = isfield(L.O, 'VicCal') && ~isempty(L.O.VicCal) || L.O.Pitch ~= 1;
    if metric
        srange = 1;
        fprintf('Metric run detected: lateral scale FIXED at 1.\n');
    else
        srange = linspace(0.5, 40, 300);
    end
    best = struct('err', inf, 's', NaN, 'sz', NaN, 'b', 0, 'sgn', 1);
    for cover_min = [0.85, 0.6, 0.4]
        for flip = 0:1
            best = scan_sb(yo, Wo, Vic, 1 - 2*flip, srange, [], cover_min, best);
        end
        if isfinite(best.s), break; end
    end
    if ~isfinite(best.s), error('compare_vic3d: global fit failed.'); end
    if metric
        best = scan_sb(yo, Wo, Vic, best.sgn, 1, ...
            linspace(best.b - 6, best.b + 6, 200), cover_min, best);
    else
        best = scan_sb(yo, Wo, Vic, best.sgn, ...
            linspace(best.s * 0.93, best.s * 1.07, 80), ...
            linspace(best.b - 6, best.b + 6, 100), cover_min, best);
    end

    % ---- report + overlay ----
    if metric
        fprintf(['\nMetric comparison (lateral scale fixed at 1): ' ...
            'depth agreement s_z = %.3f | b = %.1f mm | %s\n'], ...
            abs(best.sz), best.b, ...
            tern(best.sgn > 0, 'same direction', 'flipped'));
        if abs(abs(best.sz) - 1) > 0.10
            fprintf(['s_z should be 1.0: a residual depth discrepancy of ' ...
                '%.0f%% remains AFTER drift/tilt removal.\n'], ...
                100 * abs(abs(best.sz) - 1));
        else
            fprintf('Depth curvature agrees with Vic-3D within %.1f%%.\n', ...
                100 * abs(abs(best.sz) - 1));
        end
    else
        fprintf(['\nGlobal fit: lateral s_xy = %.3f mm | depth s_z = %.3f mm ' ...
            '| b = %.1f mm | %s\n'], best.s, abs(best.sz), best.b, ...
            tern(best.sgn > 0, 'same direction', 'flipped'));
        fprintf('Depth/lateral ratio: %.3f', abs(best.sz) / best.s);
        if abs(abs(best.sz) / best.s - 1) > 0.10
            fprintf([' <- INCONSISTENT (should be 1); with a trusted ' ...
                'calibration this points at\n   systematic effects (drift, ' ...
                'tilt) or a scale error - investigate.\n']);
        else
            fprintf(' - consistent, pitch = %.3f mm\n', best.s);
        end
    end
    fig = figure('Visible', 'off', 'Position', [40 40 950 600]); hold on;
    cols = lines(ns);
    fprintf('%-10s %9s %9s %9s %11s\n', 'state', 'rms(mm)', 'W_vic', ...
        'W_ours*s', 'drift(mm)');
    yq = best.s * (yo * best.sgn) + best.b;
    rms_state = nan(1, ns); drift_state = nan(1, ns);
    for k = 1:ns
        wv = interp1(Vic{k}.Y, Vic{k}.W, yq);
        m = isfinite(wv);
        ck = mean(wv(m)) - best.sz * mean(Wo(m, k));   % per-state offset
        wo = best.sz * Wo(:, k) + ck;
        rms_state(k)   = sqrt(mean((wo(m) - wv(m)).^2));
        drift_state(k) = -ck / best.sz;
        fprintf('%-10s %9.4f %9.4f %9.4f %11.4f\n', names{k}, ...
            rms_state(k), median(Vic{k}.W), median(wo(m)), drift_state(k));
        plot(Vic{k}.Y, Vic{k}.W, '-', 'Color', cols(k,:), ...
            'LineWidth', 1.4, 'DisplayName', [names{k} ' Vic-3D']);
        plot(yq, wo, '--', 'Color', cols(k,:), 'LineWidth', 1.2, ...
            'DisplayName', [names{k} ' LiteDIC']);
    end
    grid on; xlabel('Y, mm'); ylabel('W, mm');
    title(sprintf(['Deflection in mm: Vic-3D vs LiteDIC ' ...
        '(global fit, pitch %.3f mm)'], best.s));
    legend('Location', 'eastoutside', 'Interpreter', 'none');
    print(fig, fullfile(outDir, 'vic3d_compare.png'), '-dpng', '-r120');
    close(fig);

    % ---- how much of the reference line we actually measured ----
    % Catches a half-measured specimen; everything else is computed on the
    % overlap only.
    lo = max(min(yq), min(Vic{1}.Y));
    hi = min(max(yq), max(Vic{1}.Y));
    cover = max(0, hi - lo) / span_vic;
    fprintf(['\nCoverage: our window is Y %.1f..%.1f mm of Vic''s %.1f..%.1f ' ...
        '-> %.0f%% of the reference line\n'], min(yq), max(yq), ...
        min(Vic{1}.Y), max(Vic{1}.Y), 100 * cover);

    % ---- curvature on the window BOTH sides cover ----
    % The comparable version of the raw table above: same beam, same length.
    fprintf('\nCurvature on the matched window (detrended std, mm)\n');
    fprintf('%-10s %10s %10s %8s\n', 'state', 'vic', 'ours', 'ratio');
    curv_vic = nan(1, ns); curv_ours = nan(1, ns);
    for k = 1:ns
        inw = Vic{k}.Y >= lo & Vic{k}.Y <= hi;
        ow  = yq >= lo & yq <= hi & isfinite(Wo(:, k));
        if nnz(inw) < 5 || nnz(ow) < 5, continue; end
        curv_vic(k)  = detrended_std(Vic{k}.Y(inw), Vic{k}.W(inw));
        curv_ours(k) = detrended_std(yq(ow), best.sz * Wo(ow, k));
        fprintf('%-10s %10.4f %10.4f %8.3f\n', names{k}, curv_vic(k), ...
            curv_ours(k), curv_vic(k) / curv_ours(k));
    end

    fprintf(['\nFor metric output rerun:\n' ...
        '  run_stereo_dataset(..., ''Pitch'', %.3f)\n'], best.s);
    fprintf('Overlay saved: %s\n', fullfile(outDir, 'vic3d_compare.png'));

    % ---- BINDING verdict: speckled rows only (see header) ----
    % Vic can only measure on painted bands; bare-zone rows compare our
    % measurement against Vic's interpolation and are excluded here.
    sz_spk = NaN; rms_spk = nan(1, ns); n_spk = nnz(spk);
    if n_spk >= 10
        bs0 = struct('err', inf, 's', best.s, 'sz', NaN, 'b', best.b, ...
            'sgn', best.sgn);
        if metric
            bs0 = scan_sb(yo(spk), Wo(spk, :), Vic, best.sgn, 1, ...
                linspace(best.b - 6, best.b + 6, 200), 0.5, bs0);
        else
            bs0 = scan_sb(yo(spk), Wo(spk, :), Vic, best.sgn, ...
                linspace(best.s * 0.93, best.s * 1.07, 80), ...
                linspace(best.b - 6, best.b + 6, 100), 0.5, bs0);
        end
        if isfinite(bs0.sz)
            sz_spk = abs(bs0.sz);
            yqs = bs0.s * (yo(spk) * bs0.sgn) + bs0.b;
            Ws  = Wo(spk, :);
            for k = 1:ns
                wv = interp1(Vic{k}.Y, Vic{k}.W, yqs);
                m = isfinite(wv);
                ck = mean(wv(m)) - bs0.sz * mean(Ws(m, k));
                rms_spk(k) = sqrt(mean((bs0.sz * Ws(m, k) + ck - wv(m)).^2));
            end
        end
    end
    % COVER_MIN is part of the verdict, not a warning (DIC 3D dataset: 47 %
    % overlap passed all other checks, s_z 0.993, worst 0.6 um). A reference
    % line that overhangs the field of view fails it too - crop the export.
    % (Not to be confused with SCAN_SB's cover_min, which is how much of OUR
    % profile has to land on the reference for a candidate fit to be
    % admissible. This one is how much of the REFERENCE we measured at all.)
    REF_COVER_MIN = 0.80;
    agreed = isfinite(sz_spk) && abs(sz_spk - 1) <= 0.05 && ...
        max(rms_spk) <= 0.005;
    passed = agreed && cover >= REF_COVER_MIN;
    fprintf(['\nSPECKLE-ZONE VERDICT (rows Vic actually measures, %d rows): ' ...
        's_z = %.3f | worst rms %.4f mm | coverage %.0f%% -> %s\n'], ...
        n_spk, sz_spk, max(rms_spk), 100 * cover, tern(passed, 'PASS', 'FAIL'));
    if agreed && ~passed
        fprintf(['  the numbers agree, but only over %.0f%% of the reference ' ...
            'line (need %.0f%%):\n  the rest of the specimen was not ' ...
            'measured, and nothing above can see that.\n'], ...
            100 * cover, 100 * REF_COVER_MIN);
    end
    if ~agreed && isfinite(sz_spk)
        fprintf(['  (all-row s_z above mixes measured data with Vic''s ' ...
            'interpolated bridge over bare zones - diagnose per zone.)\n']);
    end

    % ---- programmatic verdict (identical to the printed report) ----
    out = struct();
    out.s_z         = abs(best.sz);   % depth agreement, 1.0 = perfect
    out.s_xy        = best.s;         % lateral scale / pitch, mm (1 if metric)
    out.b           = best.b;         % line-origin offset, mm
    out.orientation = best.sgn;       % +1 same direction, -1 flipped
    out.metric      = metric;         % true if VicCal or Pitch~=1
    out.rms         = rms_state;      % per-state overlay rms, mm
    out.rms_max     = max(rms_state); % worst state, mm
    out.drift       = drift_state;    % per-state whole-rig drift, mm
    out.states      = names;          % state stems (cell)
    out.s_z_speckle     = sz_spk;     % depth agreement, speckled rows only
    out.rms_speckle     = rms_spk;    % per-state rms, speckled rows, mm
    out.rms_max_speckle = max(rms_spk);
    out.n_speckle       = n_spk;      % speckled rows used
    out.cover           = cover;      % fraction of Vic's line our window spans
    out.curv_vic        = curv_vic;   % detrended std on the matched window, mm
    out.curv_ours       = curv_ours;
    out.agreed          = agreed;     % the numbers agree, coverage aside
    out.pass            = passed;     % binding verdict: agreement AND coverage
end

function s = detrended_std(y, w)
% Curvature amplitude: the spread left after removing the best straight
% line. Rigid translation and rotation of the whole specimen drop out, so
% what remains is the bending itself.
    y = y(:); w = w(:);
    Q = [ones(numel(y), 1), y];
    s = std(w - Q * (Q \ w));
end

function f = aux_file(dataDir, name)
% Locate a reference companion file (Vic export, project.xml): dataDir
% first, then its parent. The 'CalDir' layout puts the images in
% <set>\Dataset\ while the exports stay in <set>\; the flat layout keeps
% both together. Returns the dataDir path when neither exists, so callers
% that guard with isfile() keep their old behaviour.
    f = fullfile(dataDir, name);
    if isfile(f), return; end
    up = fullfile(fileparts(dataDir), name);
    if isfile(up)
        f = up;
        fprintf('compare_vic3d: using %s from the parent folder\n', name);
    end
end

function spk = speckle_rows(dataDir, L, gyk)
% Which profile rows lie on a PAINTED (speckled) band of the specimen:
% local texture (std) of the reference cam-0 frame inside the specimen
% ROI, sampled in a one-grid-step band around each kept grid row. Speckle
% paint carries several times the texture of bare metal, so a single
% relative threshold splits them cleanly. Fail-open: if the reference
% image cannot be found, every row counts as speckled.
    spk = true(size(gyk));
    stem = L.test_stems{L.sel(1)};
    dd = dir(fullfile(dataDir, [stem '_0.*']));
    if isempty(dd), return; end
    img = imread(fullfile(dataDir, dd(1).name));
    img = double(img);
    if ndims(img) == 3, img = mean(img, 3); end
    img = img / max(img(:));
    [Hi, Wi] = size(img);
    c1 = 1; c2 = Wi;
    if isfield(L.P, 'roi') && ~isempty(L.P.roi)
        c1 = max(1, round(L.P.roi(1)));
        c2 = min(Wi, round(L.P.roi(1) + L.P.roi(3)));
    end
    hb = max(4, round(L.P.grid_step / 2));
    tex = nan(size(gyk));
    for r = 1:numel(gyk)
        r1 = max(1, round(gyk(r)) - hb);
        r2 = min(Hi, round(gyk(r)) + hb);
        band = img(r1:r2, c1:c2);
        tex(r) = std(band(:));
    end
    ts = sort(tex(isfinite(tex)));
    if isempty(ts), return; end
    thr = 0.5 * ts(max(1, round(0.90 * numel(ts))));   % 0.5 x p90
    spk = tex > thr;
end

function best = scan_sb(yo, Wo, Vic, sgn, srange, bfix, cover_min, best)
% Global scan over (s_xy, b): for each candidate the depth scale s_z is
% the least-squares solution of  s_z * Wo == Wv  pooled over ALL states
% (signed - axis conventions may differ), so the error measures agreement
% of profile SHAPE and load evolution, not a single amplitude corner.
    ns = numel(Vic);
    for s = srange
        yq = s * (yo * sgn);
        if isempty(bfix)
            yv1 = Vic{1}.Y;
            bs = linspace(min(yv1) - max(yq), max(yv1) - min(yq), 150);
        else
            bs = bfix;
        end
        for b = bs
            % Per-state uniform OFFSETS are allowed on our side (whole-rig
            % drift toward the cameras rides on dZ; Vic's export is
            % drift-free), so both profiles are centered within the
            % overlap. Slope is NOT removed: in a tension test the
            % out-of-plane profile is close to linear, the slope is the signal.
            num = 0; den = 0; C = cell(1, ns); ok = true;
            for k = 1:ns
                wv = interp1(Vic{k}.Y, Vic{k}.W, yq + b);
                m = isfinite(wv);
                if nnz(m) < cover_min * numel(yq), ok = false; break; end
                woc = Wo(m, k) - mean(Wo(m, k));
                wvc = wv(m)    - mean(wv(m));
                num = num + woc.' * wvc;
                den = den + woc.' * woc;
                C{k} = {woc, wvc};
            end
            if ~ok || den < eps, continue; end
            sz = num / den;                 % signed depth scale, mm
            err = 0;
            for k = 1:ns
                e = sz * C{k}{1} - C{k}{2};
                err = err + mean(e.^2);
            end
            err = err / ns;
            if err < best.err
                best = struct('err', err, 's', s, 'sz', sz, 'b', b, 'sgn', sgn);
            end
        end
    end
end

function c = corr_safe(a, b)
    m = isfinite(a) & isfinite(b);
    if nnz(m) < 3, c = 0; return; end
    a = a(m) - mean(a(m)); b = b(m) - mean(b(m));
    d = norm(a) * norm(b);
    if d < eps, c = 0; else, c = (a.' * b) / d; end
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end

function B = parse_vic_csv(f)
% Parse the Vic-3D line export: blocks of
%   "<frame>"\n "L0";\n "Y [mm]";"W [mm]"\n  <num>;<num> ...
    txt = fileread(f);
    lines = regexp(txt, '\r?\n', 'split');
    B = struct('stem', {}, 'Y', {}, 'W', {});
    cur = [];
    for i = 1:numel(lines)
        s = strtrim(lines{i});
        if isempty(s), continue; end
        h = regexp(s, '^"(21-[\w-]+)"$', 'tokens', 'once');
        if ~isempty(h)
            cur = numel(B) + 1;
            B(cur).stem = h{1}; B(cur).Y = []; B(cur).W = [];
            continue;
        end
        if isempty(cur) || s(1) == '"', continue; end
        v = sscanf(strrep(s, ',', '.'), '%f;%f');
        if numel(v) == 2
            B(cur).Y(end+1, 1) = v(1);
            B(cur).W(end+1, 1) = v(2);
        end
    end
end
