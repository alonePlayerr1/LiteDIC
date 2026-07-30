function h = track_overlay(img, X, Y, valid, info, outfile)
% TRACK_OVERLAY  Save/show one tracking-diagnostic frame: the image with
% the tracked grid points drawn on top (green = tracked, red = lost), an
% optional ROI box and a caption. Written once per processed frame by
% RUN_SERIES_2D / RUN_STEREO_DATASET (the 'Frames' option), the sequence
% is a frame-by-frame storyboard of HOW the tracking follows the specimen,
% so a failure (grid drifting off the strip, points dying under a weld
% glow or motion blur, a wrong coarse shift) is visible by eye instead of
% being buried inside a black-box result. With info.fig it doubles as a
% LIVE view: the same overlay is drawn into a persistent visible figure
% that updates while the run progresses (RUN_SERIES_2D 'Live' option).
%
%   h = TRACK_OVERLAY(img, X, Y, valid, info, outfile)
%     img      grayscale (or RGB) image of the CURRENT frame, any scale
%     X, Y     current tracked point positions, px (any shape)
%     valid    logical, same size as X: true = point still tracked (green)
%     info     struct with optional fields:
%                .roi      [x y w h] box to draw (current specimen ROI)
%                .title    caption printed above the image
%                .lost     logical like valid: points to mark as lost
%                          (red x) at their X,Y positions
%                .lost_xy  [n x 2] explicit lost-point positions (used when
%                          a lost point has no current position, e.g. the
%                          stereo driver marks them at their reference
%                          location); overrides .lost
%                .fig      existing figure handle: draw there, keep it
%                          open and VISIBLE, drawnow (live mode). Without
%                          it a throwaway off-screen figure is used.
%     outfile  output PNG path, or '' to skip saving (live-only call)
%     h        the figure handle (reuse as info.fig on the next call)
%
% The figure is sized to the IMAGE aspect ratio (long side ~1150 pt) and
% printed at 150 dpi, so a 4k frame yields a ~2400-px-wide PNG - grid
% points stay sharp instead of being crushed into a 720p-ish thumbnail.
% Points with a non-finite position are skipped silently. Base MATLAB,
% no toolboxes.

    if nargin < 5 || isempty(info), info = struct(); end
    if nargin < 6, outfile = ''; end

    g = img;
    if ndims(g) == 3 && size(g, 3) >= 3
        g = 0.2989 * double(g(:,:,1)) + 0.5870 * double(g(:,:,2)) ...
          + 0.1140 * double(g(:,:,3));
    else
        g = double(g);
    end
    lo = min(g(:)); hi = max(g(:));
    if hi > lo, g = (g - lo) / (hi - lo); end

    % figure sized to the image aspect: long side ~1150 pt, 60 pt for title
    [Hi, Wi] = size(g);
    if Wi >= Hi, fw = 1150; fh = round(fw * Hi / Wi) + 60;
    else,        fh = 1150; fw = round((fh - 60) * Wi / Hi);
    end

    live = isfield(info, 'fig') && ~isempty(info.fig) && ishghandle(info.fig);
    if live
        f = info.fig;
        clf(f);
        figure(f);                       % bring to front, make current
    else
        f = figure('Visible', 'off', 'Position', [50 50 fw fh]);
    end

    imagesc(g); colormap(gray); axis image; hold on;

    ok = valid(:) & isfinite(X(:)) & isfinite(Y(:));
    plot(X(ok), Y(ok), '.', 'Color', [0 0.9 0], 'MarkerSize', 4);

    if isfield(info, 'lost_xy') && ~isempty(info.lost_xy)
        lx = info.lost_xy(:, 1); ly = info.lost_xy(:, 2);
        plot(lx, ly, 'x', 'Color', [1 0.1 0.1], 'MarkerSize', 5, 'LineWidth', 0.8);
    elseif isfield(info, 'lost') && any(info.lost(:))
        lm = info.lost(:) & isfinite(X(:)) & isfinite(Y(:));
        plot(X(lm), Y(lm), 'x', 'Color', [1 0.1 0.1], 'MarkerSize', 5, 'LineWidth', 0.8);
    end

    if isfield(info, 'roi') && numel(info.roi) == 4
        r = info.roi;
        plot(r(1) + [0 r(3) r(3) 0 0], r(2) + [0 0 r(4) r(4) 0], ...
            '-', 'Color', [1 1 0.2], 'LineWidth', 1.0);
    end

    if isfield(info, 'title') && ~isempty(info.title)
        title(info.title, 'Interpreter', 'none', 'FontSize', 9);
    end
    set(gca, 'XTick', [], 'YTick', []);
    hold off;

    if ~isempty(outfile)
        print(f, outfile, '-dpng', '-r150');
    end
    if live
        drawnow;                         % update the on-screen view
    else
        close(f); f = gobjects(0);
    end
    if nargout > 0, h = f; end
end
