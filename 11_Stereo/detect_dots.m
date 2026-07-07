function D = detect_dots(img, params)
% DETECT_DOTS  Toolbox-free detection of a black-on-white dot calibration
% target with "donut" fiducial dots (bright-centre markers).
%
%   D = DETECT_DOTS(img) takes a grayscale (or RGB) image of a regular grid
%   of black dots on a light background and returns the subpixel dot
%   centres and which dots are donut fiducials. No toolboxes are used.
%
%   Output struct D:
%     D.x, D.y     subpixel dot centres, px (column vectors, x = column)
%     D.is_donut   logical column, true for the fiducial (donut) dots
%     D.donut_xy   [k x 2] donut centres [x y], k is usually 3
%     D.n          number of detected dots
%     D.diam       median dot diameter, px
%     D.spacing    median nearest-neighbour spacing, px
%
%   params (optional struct, all fields have defaults):
%     .bg_win    background window, px, must exceed the dot diameter (61)
%     .rel_thr   relative-darkness threshold for a dot pixel, 0..1 (0.30)
%     .donut_thr centre-brightness ratio above which a dot is a donut (0.30)
%     .min_fill  minimum area / bounding-box-area fill ratio (0.45)
%     .area_lo   low  area cut as a fraction of the median dot area (0.40)
%     .area_hi   high area cut as a fraction of the median dot area (2.50)
%     .border    reject blobs within this many px of the image edge (2)
%     .board_win window for the board (paper) mask, px (81)
%     .board_frac board-mask level as a fraction of the 90th-pct brightness
%                (0.32); masks out the dark surround and its boundary halo
%
%   Method: the local paper brightness is estimated with a separable
%   moving maximum over a window larger than one dot (the same device
%   dic_engine.m uses for its "hole" mask), which removes the strong
%   vignetting of these frames. A pixel is part of a dot when it is darker
%   than its local paper level by rel_thr. Dark pixels are grouped into
%   connected components by iterative label-min propagation (4-connected),
%   filtered to dot-sized round blobs that do not touch the border, and
%   reduced to darkness-weighted subpixel centres. A donut is told from a
%   solid dot by the brightness at its centre: a solid dot is darkest at
%   the centre, a donut is paper-bright there.
%
%   See also DIC_ENGINE (movmax hole mask), INDEX_GRID.

    if nargin < 2, params = struct(); end
    P = defaults(params);

    g = to_gray(img);
    [H, W] = size(g);

    % --- background flatten: local paper brightness, then relative darkness
    bg  = movmax(movmax(g, P.bg_win, 1), P.bg_win, 2);
    rel = (bg - g) ./ max(bg, 1);
    % --- board mask: keep only where the neighbourhood is bright paper on
    % average. This drops the dark mount/surround AND the ~bg_win-wide halo
    % band along the board boundary, where the moving-max background is
    % bright but the pixels themselves are dark surround. That band is one
    % large connected region that would otherwise dominate - and drastically
    % slow - component labeling on views with a big dark border.
    paper = prctile_(g, 90);
    lm    = movmean(movmean(g, P.board_win, 1), P.board_win, 2);
    board = lm > P.board_frac * paper;
    bw    = (rel > P.rel_thr) & board;

    % --- connected components (4-connected) ---
    [L, n] = cc_label(bw);
    if n == 0, D = empty_result(); return; end
    A = accumarray(L(L > 0), 1);

    % --- keep dot-sized, round, non-border blobs ---
    [x0, x1, y0, y1] = blob_bbox(L, n);
    bbarea = (x1 - x0 + 1) .* (y1 - y0 + 1);
    med = median(A(A >= 200));
    if ~isfinite(med) || med <= 0, med = median(A); end
    touches = x0 <= P.border | y0 <= P.border | x1 >= W - P.border | y1 >= H - P.border;
    fill = A ./ max(bbarea, 1);
    isdot = A > P.area_lo * med & A < P.area_hi * med & ~touches & fill > P.min_fill;
    dots = find(isdot);
    if isempty(dots), D = empty_result(); return; end

    % --- darkness-weighted subpixel centres ---
    dk = max(bg - g, 0);
    nd = numel(dots);
    cx = zeros(nd, 1); cy = zeros(nd, 1);
    for i = 1:nd
        m = (L == dots(i));
        w = dk(m); [rr, cc] = find(m);
        sw = sum(w);
        cx(i) = sum(cc .* w) / sw;
        cy(i) = sum(rr .* w) / sw;
    end

    % --- donut fiducials by centre brightness (robust to ring breaks) ---
    rad = max(4, round(0.15 * 2 * sqrt(med / pi)));   % ~15% of dot diameter
    cscore = zeros(nd, 1);
    for i = 1:nd
        xi = round(cx(i)); yi = round(cy(i));
        ys = max(1, yi - rad):min(H, yi + rad);
        xs = max(1, xi - rad):min(W, xi + rad);
        patch = g(ys, xs);
        cscore(i) = mean(patch(:)) / max(bg(yi, xi), 1);
    end
    % The target has exactly three donut fiducials. Keep the three strongest
    % bright-centre dots above the threshold; capping at three rejects the
    % occasional solid dot with a bright blemish at its centre, while a view
    % that yields fewer than three fails safe (INDEX_GRID drops it).
    is_donut = false(nd, 1);
    [~, order] = sort(cscore, 'descend');
    cand = order(cscore(order) > P.donut_thr);
    is_donut(cand(1:min(3, numel(cand)))) = true;

    D.x = cx; D.y = cy;
    D.is_donut = is_donut;
    D.donut_xy = [cx(is_donut), cy(is_donut)];
    D.n = nd;
    D.diam = 2 * sqrt(median(A(dots)) / pi);
    D.spacing = nn_spacing(cx, cy);
    D.cscore = cscore;
end

% ===================== helpers (all base MATLAB) =====================

function P = defaults(params)
    d = struct('bg_win', 61, 'rel_thr', 0.30, 'donut_thr', 0.30, ...
        'min_fill', 0.45, 'area_lo', 0.40, 'area_hi', 2.50, 'border', 2, ...
        'board_win', 81, 'board_frac', 0.32);
    fn = fieldnames(d);
    for k = 1:numel(fn)
        if isfield(params, fn{k}) && ~isempty(params.(fn{k}))
            d.(fn{k}) = params.(fn{k});
        end
    end
    P = d;
end

function g = to_gray(img)
    if ndims(img) == 3 && size(img, 3) >= 3
        img = double(img);
        g = 0.2989 * img(:,:,1) + 0.5870 * img(:,:,2) + 0.1140 * img(:,:,3);
    else
        g = double(img);
    end
end

function v = prctile_(x, q)
    x = sort(x(isfinite(x))); n = numel(x);
    if n == 0, v = 0; return; end
    v = x(min(n, max(1, round(q / 100 * n))));
end

function [L, n] = cc_label(bw)
% 4-connected components via union-find over the dark pixels. This is
% O(#dark pixels) and independent of blob geometry: a large stray region
% (e.g. a residual boundary band) costs no more than its pixel count, so
% labeling time never blows up on strongly tilted views the way an
% iterative flood/label-propagation would.
    [H, W] = size(bw);
    lin = find(bw); M = numel(lin);
    L = zeros(H, W);
    if M == 0, n = 0; return; end
    nid = zeros(H*W, 1); nid(lin) = 1:M;          % node id per dark pixel
    r = mod(lin - 1, H) + 1;                       % row (column-major)
    c = floor((lin - 1) / H) + 1;                  % col
    % edges to the right (+H) and down (+1) neighbours that are also dark
    src = nid(lin(c < W)); nb = nid(lin(c < W) + H);
    ER = [src(nb > 0), nb(nb > 0)];
    src = nid(lin(r < H)); nb = nid(lin(r < H) + 1);
    ED = [src(nb > 0), nb(nb > 0)];
    E = [ER; ED];
    parent = (1:M)';
    for e = 1:size(E, 1)
        a = E(e,1); b = E(e,2);
        while parent(a) ~= a, parent(a) = parent(parent(a)); a = parent(a); end
        while parent(b) ~= b, parent(b) = parent(parent(b)); b = parent(b); end
        if a ~= b, parent(a) = b; end
    end
    root = zeros(M, 1);
    for k = 1:M
        a = k; while parent(a) ~= a, a = parent(a); end
        root(k) = a;
    end
    [~, ~, g] = unique(root); n = max(g);
    L(lin) = g;
end

function [x0, x1, y0, y1] = blob_bbox(L, n)
    x0 = zeros(n,1); x1 = x0; y0 = x0; y1 = x0;
    [rr, cc] = find(L > 0);
    lab = L(L > 0);
    x0 = accumarray(lab, cc, [n 1], @min);
    x1 = accumarray(lab, cc, [n 1], @max);
    y0 = accumarray(lab, rr, [n 1], @min);
    y1 = accumarray(lab, rr, [n 1], @max);
end

function d = nn_spacing(cx, cy)
    N = numel(cx);
    if N < 2, d = NaN; return; end
    dmin = inf(N, 1);
    for i = 1:N
        dd = hypot(cx - cx(i), cy - cy(i)); dd(i) = inf;
        dmin(i) = min(dd);
    end
    d = median(dmin);
end

function D = empty_result()
    D = struct('x', [], 'y', [], 'is_donut', logical([]), 'donut_xy', [], ...
        'n', 0, 'diam', NaN, 'spacing', NaN, 'cscore', []);
end
