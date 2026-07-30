function export_profiles(dataDir)
% EXPORT_PROFILES  Dump the stereo deflection profiles W(Y) to a plain CSV.
%
%   export_profiles('C:\...\<dataset>')  reads <dataDir>\litedic_out\
%   stereo_results.mat (produced by RUN_STEREO_DATASET) and writes
%   <dataDir>\litedic_out\profiles_wy.csv with one row per (state, grid
%   row): state stem, along-beam coordinate Y (units of the run), and the
%   per-row median normal deflection W. The same plane frame (c0, e1, nrm)
%   the driver saved is reused, so the numbers are identical to the ones
%   compare_vic3d consumes - this is the raw-data export for external
%   plotting, digging into a validation mismatch, or archiving alongside a
%   paper figure. Base MATLAB, no toolboxes.
    outDir = fullfile(dataDir, 'litedic_out');
    L = load(fullfile(outDir, 'stereo_results.mat'));
    [nY, nX, ~] = size(L.S.X_ref_mm);
    P3 = reshape(L.S.X_ref_mm, [], 3);
    Sc = reshape((P3 - L.c0) * L.e1, nY, nX);   % along-beam coordinate

    fout = fullfile(outDir, 'profiles_wy.csv');
    fid = fopen(fout, 'w');
    fprintf(fid, 'state,y,W\n');
    for k = 2:numel(L.sel)
        d3 = reshape(L.S.frames(k).d3D_mm, [], 3);
        Wn = reshape(d3 * L.nrm, nY, nX);
        for r = 1:nY
            yr = median(Sc(r, isfinite(Wn(r, :))), 'omitnan');
            wr = median(Wn(r, :), 'omitnan');
            if isfinite(yr) && isfinite(wr)
                fprintf(fid, '%s,%.6f,%.8f\n', ...
                    L.test_stems{L.sel(k)}, yr, wr);
            end
        end
    end
    fclose(fid);
    fprintf('Saved %s\n', fout);
end
