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
%             .clim   [lo hi] color limits (default: the tighter of
%                     median +/- 3 MAD-sigma and the 1..99 percentiles - see
%                     rlim() below; same rule as RUN_SERIES_2D and
%                     RUN_STEREO_DATASET)
%             .levels number of discrete colour bands (default 30, as Vic-3D; 0 = a
%                     continuous map). Banded is how every DIC package
%                     draws a field, and it is what makes a LiteDIC map and
%                     a Vic-3D map of the same specimen READ alike: on a
%                     continuous map the eye judges the gradient, and two
%                     equal fields with different limits look different.
%             .cmap   colormap name or N x 3 matrix (default 'dic_rainbow':
%                     magenta-blue-cyan-green-yellow-red in Vic-3D's order,
%                     resampled for even perceptual steps - see DIC_RAINBOW;
%                     plain jet/hsv show false bands in cyan and yellow)
%             .alpha  opacity of the colour over the frame (default 0.62,
%                     about what Vic-3D uses: the speckle stays visible
%                     through the map, so the eye can tell field from photo)
%             .smooth display-only Gaussian smoothing of a GRIDDED field,
%                     in grid nodes (default 1; 0 = off). Node noise of a few
%                     microns jitters across 0.1 mm band edges and makes
%                     them ragged; one node of smoothing removes the jitter.
%                     The numbers in the .mat are untouched.
%
%   When X, Y, V are matrices of the same size (a tracked grid) the field is
%   drawn as filled contour bands - smooth edges, like every DIC package's
%   map. Vectors fall back to coloured points.
%             .ref    'raw' (default) or 'plane': subtract the least-squares
%                     plane V = a + b X + c Y fitted through the drawn nodes
%                     before colouring, and put the plane's parameters in the
%                     title. This is how Vic-3D draws its deflection map ("W
%                     from best-fit plane"): the same bowl, referenced to the
%                     reference frame, spans -2.1..-0.3 mm on our raw map and
%                     -0.8..+0.9 mm on theirs, and side by side the two colour
%                     bars looked like a disagreement. Same convention, same
%                     scale. The offset and tilt of the plane are the only
%                     things removed; the numbers in the .mat are untouched.
%             .hide_outliers  true (default): nodes further outside the
%                     colour limits than 10 % of the range are not drawn.
%                     They are triangulation failures at the mask edge -
%                     on the sandwich series one sat at -72 mm on a 10 mm
%                     bowl - not measurement, and drawn saturated at the
%                     rim they made a clean field look like rubbish. The
%                     count goes into the title so nothing is hidden
%                     silently; the numbers in the .mat are untouched.
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
    levels = 30;
    if isfield(opts, 'levels') && ~isempty(opts.levels), levels = opts.levels; end
    cmap = 'dic_rainbow';
    if isfield(opts, 'cmap') && ~isempty(opts.cmap), cmap = opts.cmap; end
    alpha_ = 0.62;
    if isfield(opts, 'alpha') && ~isempty(opts.alpha), alpha_ = opts.alpha; end
    hide = true;
    if isfield(opts, 'hide_outliers') && ~isempty(opts.hide_outliers)
        hide = logical(opts.hide_outliers);
    end

    user_clim = isfield(opts, 'clim') && numel(opts.clim) == 2;
    to_plane = isfield(opts, 'ref') && strcmpi(char(opts.ref), 'plane');
    % A caller's clim refers to the values that get COLOURED. With ref =
    % 'plane' those are plane-referenced, so the outlier rule below must not
    % compare the raw values against it (a -2.1..-0.1 mm bowl against a
    % -0.7..0.9 mm scale would hide almost every node): raw outliers are
    % then judged by the robust rule and the caller's clim is applied after
    % the plane is subtracted.
    if user_clim && ~to_plane
        c = opts.clim;
    else
        c = rlim(V(m));
    end
    nout = 0;
    if hide && c(2) > c(1)
        pad = 0.1 * (c(2) - c(1));
        out = m & (V(:) < c(1) - pad | V(:) > c(2) + pad);
        nout = nnz(out);
        m = m & ~out;
    end
    plane = [];
    if to_plane && nnz(m) >= 3
        A = [ones(nnz(m), 1), X(m), Y(m)];
        plane = A \ V(m);                       % least squares, drawn nodes only
        V = V - (plane(1) + plane(2) * X + plane(3) * Y);
        if user_clim, c = opts.clim; else, c = rlim(V(m)); end
    end

    sm = 1;
    if isfield(opts, 'smooth') && ~isempty(opts.smooth), sm = opts.smooth; end
    gridded = isequal(size(X), size(Y), size(V)) && ismatrix(V) && all(size(V) >= 3);
    if gridded && c(2) > c(1) && levels > 0
        Vm = reshape(m, size(V));
        Vg = V;  Vg(~Vm) = NaN;
        if sm > 0, Vg = smooth_nan(Vg, sm); end
        % contourf needs finite, well-behaved coordinates everywhere. Nodes
        % that lost tracking have NaN positions; a few tracked nodes have a
        % position spike (a subset that latched onto the wrong speckle) and
        % a plausible value - drawn at the wrong place they stretch a cell
        % across the whole image as a coloured streak. Both get the local
        % (5 x 5 median) position; a spike further than one grid pitch from
        % it also loses its value, so it is simply not drawn.
        [Xg, Yg, spike] = tidy_grid(X, Y);
        Vg(spike) = NaN;
        edges = linspace(c(1), c(2), levels + 1);
        keep = isfinite(Vg);                 % MATLAB's max/min drop NaN: keep the holes
        Vg = min(max(Vg, c(1)), c(2));  Vg(~keep) = NaN;
        [~, hc] = contourf(Xg, Yg, Vg, edges, 'LineStyle', 'none');
        try, hc.FaceAlpha = alpha_; catch, end   % translucent, the speckle shows through
    else
        scatter(X(m), Y(m), ms, V(m), 'filled');
    end
    cb = colorbar;                                       %#ok<NASGU>
    if c(2) > c(1), clim(c); end
    if ischar(cmap) || isstring(cmap)
        if levels > 0, map = feval(char(cmap), levels); else, map = feval(char(cmap)); end
    else
        map = cmap;
    end
    colormap(gca, map);

    ttl = '';
    if isfield(opts, 'title') && ~isempty(opts.title), ttl = opts.title; end
    if nout > 0
        ttl = sprintf('%s  (%d outlier node%s not drawn)', ttl, nout, tern(nout > 1, 's', ''));
    end
    if ~isempty(plane)
        ttl = sprintf('%s  [minus best-fit plane: %.3f %+.2e X %+.2e Y]', ttl, plane(1), plane(2), plane(3));
    end
    if ~isempty(ttl)
        title(ttl, 'Interpreter', 'tex', 'FontSize', 10);
    end
    set(gca, 'XTick', [], 'YTick', []);
    hold off;
end

function c = rlim(A)
% Robust color limits: the TIGHTER of two rules in each direction - median
% +/- 3 robust-sigma (MAD) and the 1st..99th percentiles of the data.
%
% MAD alone is right for a strain field (a narrow peak with a heavy tail of
% ROI-edge nodes) and wrong for a smooth deflection bowl: there the values
% spread almost uniformly over the range, the MAD is ~0.3 of it, and +/- 3
% sigma reaches nearly a full range PAST the data on both sides. On the
% 122-state sandwich-panel series the last state's field lay in [-10.9, 0.3]
% mm and was drawn on a [-17.4, +6.4] mm colourbar: half the colours spent on
% values that do not exist, the bowl a pale wash next to Vic-3D's picture of
% the same panel. The percentile pins the limits to where the data is; the
% MAD still cuts a heavy tail a percentile alone would let through.
% Toolbox-free; degenerate input falls back sensibly.
    a = A(isfinite(A));
    if numel(a) < 5, c = [0 1]; return; end
    m = median(a);
    d = 1.4826 * median(abs(a - m));       % robust sigma from MAD
    if d <= 0, d = std(a); end             % all-equal-ish -> fall back to std
    if d <= 0, c = [m - 1, m + 1]; return; end
    c = [m - 3*d, m + 3*d];
    s = sort(a);  n = numel(s);            % 1st / 99th percentile, toolbox-free
    p = s(max(1, min(n, round([0.01 0.99] * n))));
    c = [max(c(1), p(1)), min(c(2), p(2))];   % the tighter rule in each direction
    if ~(c(2) > c(1)), c = [m - 3*d, m + 3*d]; end   % near-flat data: MAD alone
end

function out = tern(c, a, b)
    if c, out = a; else, out = b; end
end

function G = smooth_nan(V, s)
% Normalised Gaussian smoothing that ignores NaN (display only). Toolbox-free.
    r = ceil(3 * s);  k = exp(-((-r:r).^2) / (2 * s^2));  k = k / sum(k);
    ok = isfinite(V);
    Z = V;  Z(~ok) = 0;
    num = conv2(k, k, Z, 'same');
    den = conv2(k, k, double(ok), 'same');
    G = num ./ max(den, 1e-9);
    G(~ok) = NaN;
end

function [X, Y, spike] = tidy_grid(X, Y)
% Finite, spike-free node positions for contourf (display only).
    Xm = local_median(X);  Ym = local_median(Y);
    % grid pitch in px from the median node spacing
    dx = median(abs(diff(Xm, 1, 2)), 'all', 'omitnan');
    dy = median(abs(diff(Ym, 1, 1)), 'all', 'omitnan');
    pitch = max([dx, dy, 1]);
    spike = isfinite(X) & isfinite(Y) & (hypot(X - Xm, Y - Ym) > pitch);
    rep = spike | ~isfinite(X) | ~isfinite(Y);
    X(rep) = Xm(rep);  Y(rep) = Ym(rep);
    % anything still NaN (a 5 x 5 hole with no finite neighbour): row / column median
    X = fill_lines(X, 2);  Y = fill_lines(Y, 1);
end

function M = local_median(A)
% 5 x 5 median ignoring NaN, toolbox-free (25 shifted copies).
    [nr, nc] = size(A);
    S = nan(nr, nc, 25);  q = 0;
    for di = -2:2
        for dj = -2:2
            q = q + 1;
            B = nan(nr, nc);
            r0 = max(1, 1 - di):min(nr, nr - di);
            c0 = max(1, 1 - dj):min(nc, nc - dj);
            B(r0, c0) = A(r0 + di, c0 + dj);
            S(:, :, q) = B;
        end
    end
    M = median(S, 3, 'omitnan');
end

function A = fill_lines(A, dim)
% Replace remaining NaN coordinates by the median of their row (dim = 1) or
% column (dim = 2).
    bad = ~isfinite(A);
    if ~any(bad(:)), return; end
    med = median(A, dim, 'omitnan');
    if dim == 1, M = repmat(med, size(A, 1), 1); else, M = repmat(med, 1, size(A, 2)); end
    M(~isfinite(M)) = 0;
    A(bad) = M(bad);
end
