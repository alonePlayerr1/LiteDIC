function probe_caltarget(dataDir, stem)
% PROBE_CALTARGET  Diagnostic: can detect_dots + index_grid handle a given
% calibration-target image pair? Tries both polarities (the DICe target is
% black-on-white; Correlated Solutions / Vic-3D targets are white-on-black,
% which is handled by inverting the image) and reports what it finds.
%
%   probe_caltarget('C:\...\DIC 3D dataset', '21-0293')
%
% Prints, per camera and polarity: detected dots, median spacing, donut
% count, and whether a consistent grid was indexed. Saves overlay PNGs
% (detected centres + donuts + grid axes) to <dataDir>\litedic_out\.
% Diagnostic only - modifies nothing.

    if nargin < 2, stem = '21-0293'; end
    outDir = fullfile(dataDir, 'litedic_out');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    for cam = 0:1
        f = fullfile(dataDir, sprintf('%s_%d.tif', stem, cam));
        if ~exist(f, 'file'), error('Not found: %s', f); end
        img = double(imread(f));
        img = img / max(img(:));
        if ndims(img) == 3, img = mean(img, 3); end
        fprintf('\n=== %s_%d.tif  (%d x %d px) ===\n', stem, cam, ...
            size(img, 2), size(img, 1));

        for inv = [false true]
            if inv, g = 1 - img; tag = 'inverted (white dots)';
            else,   g = img;     tag = 'as-is (black dots)';
            end
            try
                D = detect_dots(g);
            catch ME
                fprintf('  %-22s detect_dots ERROR: %s\n', tag, ME.message);
                continue;
            end
            gok = false; gn = 0;
            if D.n >= 12 && size(D.donut_xy, 1) == 3
                try
                    G = index_grid(D);
                    gok = G.ok; gn = G.n;
                catch
                    % no usable grid: gok stays false, reported as "no"
                end
            end
            fprintf('  %-22s dots %4d | spacing %6.1f px | diam %5.1f px | donuts %d | grid %s (%d indexed)\n', ...
                tag, D.n, D.spacing, D.diam, size(D.donut_xy, 1), ...
                tern(gok, 'OK', 'no'), gn);

            % overlay for eyeballing
            fig = figure('Visible', 'off');
            imagesc(g); colormap(gray); axis image; hold on;
            plot(D.x, D.y, 'g+', 'MarkerSize', 6);
            if ~isempty(D.donut_xy)
                plot(D.donut_xy(:,1), D.donut_xy(:,2), 'ro', ...
                    'MarkerSize', 12, 'LineWidth', 1.5);
            end
            title(sprintf('%s_%d %s: %d dots, %d donuts', stem, cam, tag, ...
                D.n, size(D.donut_xy, 1)), 'Interpreter', 'none');
            print(fig, fullfile(outDir, sprintf('probe_%s_%d_%s.png', ...
                stem, cam, tern(inv, 'inv', 'raw'))), '-dpng', '-r100');
            close(fig);
        end
    end
    fprintf('\nOverlays saved to %s\n', outDir);
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end
