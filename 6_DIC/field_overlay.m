function field_overlay(img, X, Y, V, opts)
% FIELD_OVERLAY  Draw a DIC field ON the specimen itself: the grayscale
% frame as an RGB background (so the axes colormap stays free for the
% data) with the field values as filled colored points at their CURRENT
% (tracked, deformed) positions.
%
% Why: the classic imagesc panels draw the Lagrangian fields on the
% frame-1 grid, which is always an axis-aligned rectangle - a bent beam
% still shows up as a straight horizontal band. Here the same values are
% plotted where the material points actually ARE in the displayed frame,
% so the color map follows the real deformed shape of the object.
%
%   FIELD_OVERLAY(img, X, Y, V, opts) draws into the CURRENT axes
%   (caller owns the figure/subplot/print):
%     img   grayscale (or RGB) frame, any scale
%     X, Y  current point positions, px (any shape)
%     V     field value per point (same size); NaN points are skipped
%     opts  optional struct:
%             .title  axes title
%             .msize  scatter marker area, pt^2 (default 6)
%             .clim   [lo hi] color limits (default: robust median +/- 3
%                     MAD-sigma, so a few outliers don't crush the scale)
%
% Base MATLAB, no toolboxes.
    if nargin < 5 || isempty(opts), opts = struct(); end

    g = img;
    if ndims(g) == 3 && size(g, 3) >= 3
        g = 0.2989 * double(g(:,:,1)) + 0.5870 * double(g(:,:,2)) ...
          + 0.1140 * double(g(:,:,3));
    else
        g = double(g);
    end
    lo = min(g(:)); hi = max(g(:));
    if hi > lo, g = (g - lo) / (hi - lo); end

    image(repmat(g, [1 1 3])); axis image; hold on;

    m = isfinite(V(:)) & isfinite(X(:)) & isfinite(Y(:));
    ms = 6;
    if isfield(opts, 'msize') && ~isempty(opts.msize), ms = opts.msize; end
    scatter(X(m), Y(m), ms, V(m), 'filled');
    cb = colorbar;                                       %#ok<NASGU>

    if isfield(opts, 'clim') && numel(opts.clim) == 2
        clim(opts.clim);
    else
        a = V(m);
        if numel(a) >= 5
            md = median(a);
            d = 1.4826 * median(abs(a - md));
            if d <= 0, d = std(a); end
            if d > 0, clim([md - 3*d, md + 3*d]); end
        end
    end
    if isfield(opts, 'title') && ~isempty(opts.title)
        title(opts.title, 'Interpreter', 'tex', 'FontSize', 10);
    end
    set(gca, 'XTick', [], 'YTick', []);
    hold off;
end
