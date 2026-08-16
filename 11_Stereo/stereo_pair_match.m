function D = stereo_pair_match(img0, img1, P)
% STEREO_PAIR_MATCH  Cross-camera subset matching with the EXISTING 2D
% IC-GN engine: a global integer pre-alignment absorbs the large stereo
% disparity offset, then DIC_ENGINE (untouched) resolves the residual
% field with its reliability-guided propagation.
%
%   D = STEREO_PAIR_MATCH(img0, img1, P) matches the synchronized frames
%   img0 (camera 0, used as the DIC reference) and img1 (camera 1), both
%   grayscale double in [0..1]. The total left->right correspondence of a
%   camera-0 point (x, y) is
%       (x, y)  ->  (x + D.U(iy,ix), y + D.V(iy,ix))
%   on the grid D.gx / D.gy.
%
%   Why the pre-alignment: DIC_ENGINE initializes from ONE integer search
%   and propagates; a stereo pair has a large (tens of px) nearly-uniform
%   disparity offset plus a smooth few-px variation. Shifting img1 by the
%   global integer offset first (found by FFT cross-correlation, base
%   MATLAB) keeps the residual field small and smooth - exactly the regime
%   the engine's RG propagation is built for - and keeps the engine's
%   margin (subset/2 + search_range) small so border subsets survive.
%   The shift is a pure integer index move (no interpolation), so engine
%   numerics are unaffected; it is added back to the output verbatim.
%
%   P: passed to DIC_ENGINE (subset_size, grid_step, search_range,
%      black_thr, zncc_min, max_iter, conv_tol), plus optional
%      P.max_global (clip for the global shift sanity check, default 300 px),
%      P.shift_prior ([dx dy] expected pre-alignment - see below) and
%      P.shift_win ([wx wy] half-window around that prior, default [40 4]).
%
%   Output struct D:
%     D.U, D.V    total disparity components, px (NaN where invalid)
%     D.gx, D.gy  camera-0 grid coordinates (same convention as DIC_ENGINE)
%     D.Q, D.NIT  engine quality / iteration fields
%     D.shift     the global integer [dx dy] that was pre-applied
%     D.peakq     normalized height of the correlation peak (diagnostic)
%     D.shift_prior  the prior used ([] when the peak was picked globally)
%     D.shift_ok  false when no credible peak existed inside the window and
%                 the prior itself was used as the pre-alignment
%
%   Why the prior: the unconstrained argmax of the correlation is not safe.
%   Measured on the sandwich sets, 12 of 129 frames of series 11 (the
%   specimen's failure event, 11-0087..11-0091 and 11-0097..11-0103) put the
%   global maximum at offsets like (-300, 1068) px - with a normalized peak
%   height INSIDE the range the healthy frames occupy (0.56-0.74 vs
%   0.61-0.78), so no peak-quality gate can separate them. One of those
%   exceeds max_global and aborts the whole series from here; the rest pass
%   silently and give wrong states. The cross-camera offset is a RIG
%   quantity, though. Measured over all 230 frames of the three sandwich
%   sets, between consecutive frames it moves at most 25 px in dx and 1 px
%   in dy, so [prior +/- [40 4]] contains every legitimate motion with room
%   to spare while every one of the 12 failures misses it by 234 px or more
%   in dy alone.
%
%   Constrained search alone is not enough: if the true maximum is outside
%   the window, the argmax lands on the window BOUNDARY, also wrong (seen: dx pinned at
%   -151 px, one state left with a 56 px epipolar residual). So a boundary
%   hit is treated as "no credible peak here" and the prior is used verbatim
%   instead - safe, because the pre-alignment is only a warm start that
%   DIC_ENGINE resolves the remainder of, and it is added back to the output
%   unchanged. D.shift_ok reports which of the two happened.
%
%   With no prior the whole plane is searched, exactly as before.
%
%   See also DIC_ENGINE, STEREO_DIC, STEREO_GEOMETRY.

    if nargin < 3, error('stereo_pair_match:args', 'img0, img1, P required.'); end
    maxg = 300;
    if isfield(P, 'max_global') && ~isempty(P.max_global), maxg = P.max_global; end
    prior = [];
    if isfield(P, 'shift_prior') && ~isempty(P.shift_prior), prior = P.shift_prior; end
    win = [40 4];
    if isfield(P, 'shift_win') && ~isempty(P.shift_win), win = P.shift_win; end
    % dic_engine gets only its own fields
    P = rmfield_safe(P, {'max_global', 'shift_prior', 'shift_win'});

    % ---- 1) global integer offset by FFT cross-correlation ----
    [dx, dy, peakq, shift_ok] = global_shift(img0, img1, prior, win);
    if abs(dx) > maxg || abs(dy) > maxg
        error('stereo_pair_match:shift', ...
            'global shift (%d, %d) exceeds max_global = %d px.', dx, dy, maxg);
    end

    % ---- 2) shift img1 so the pair nearly overlaps ----
    % img1 feature of the camera-0 point (x,y) sits at (x+dx, y+dy);
    % build shifted(y, x) = img1(y+dy, x+dx). Uncovered borders become 0
    % (black), which the engine masks/rejects on its own.
    [H, W] = size(img1);
    shifted = zeros(H, W);
    ys = max(1, 1-dy):min(H, H-dy);      % rows where y+dy is inside img1
    xs = max(1, 1-dx):min(W, W-dx);
    shifted(ys, xs) = img1(ys+dy, xs+dx);

    % ---- 3) residual field with the untouched 2D engine ----
    R = dic_engine(img0, shifted, P);

    D.U = R.U + dx;  D.V = R.V + dy;
    D.gx = R.gx;     D.gy = R.gy;
    D.Q = R.Q;       D.NIT = R.NIT;
    D.nseeds = R.nseeds;   % how many fronts it took to fill this field;
                           % a jump between states means the field broke up
    D.nfail  = R.nfail;         % seed attempts that found nothing
    D.seed_capped = R.seed_capped;  % true = the budget ran out, not the grid
    D.unreached = R.unreached;          % speckle fraction never propagated to
    D.unreached_max = R.unreached_max;  % ... its largest connected piece
    D.shift = [dx dy];
    D.peakq = peakq;
    D.shift_prior = prior;
    D.shift_ok = shift_ok;
end

% ---------------------------------------------------------------------
function [dx, dy, q, ok] = global_shift(a, b, prior, win)
% Integer offset (dx, dy) with b(y+dy, x+dx) ~ a(y, x), from the peak of
% the circular FFT cross-correlation of the zero-mean, Hann-tapered
% images. c = ifft2(fft2(b) .* conj(fft2(a))) gives
%   c(1+dy, 1+dx) = sum_x b(x + [dx dy]) .* a(x)   (indices wrap),
% so the peak location IS the offset of b relative to a.
%
% With PRIOR = [dx0 dy0] the maximum is taken over the candidate offsets
% [dx0 +/- win(1)] x [dy0 +/- win(2)] instead of the whole plane: the offset
% is a rig quantity and cannot leave that neighbourhood physically, while
% the unconstrained argmax demonstrably does (see the header). A maximum ON
% the window boundary means the real one is outside it, i.e. this frame has
% no credible global peak; then the prior is returned with ok = false rather
% than the boundary point. Without a prior the whole plane is searched, as
% before, and ok is true.
    if nargin < 3, prior = []; end
    if nargin < 4 || isempty(win), win = [40 4]; end
    ok = true;
    [H, W] = size(a);
    wy = hann_win(H); wx = hann_win(W);
    Wn = wy * wx.';                       % separable taper, kills wrap edges
    az = (a - mean(a(:))) .* Wn;
    bz = (b - mean(b(:))) .* Wn;
    c = real(ifft2(fft2(bz) .* conj(fft2(az))));
    if isempty(prior)
        [q, k] = max(c(:));
        [iy, ix] = ind2sub([H, W], k);
        dy = iy - 1; dx = ix - 1;
        if dy > H/2, dy = dy - H; end     % unwrap to signed offsets
        if dx > W/2, dx = dx - W; end
    else
        % candidate offsets around the prior, mapped onto the wrapped
        % correlation indices (offset d lives at index d+1, modulo the size)
        ox = round(prior(1)) + (-abs(win(1)) : abs(win(1)));
        oy = round(prior(2)) + (-abs(win(2)) : abs(win(2)));
        ox = ox(abs(ox) <= floor(W/2));   % offsets beyond +/-N/2 are aliases
        oy = oy(abs(oy) <= floor(H/2));
        if isempty(ox) || isempty(oy)
            error('stereo_pair_match:win', ...
                'shift_prior/shift_win select no valid offset.');
        end
        sub = c(mod(oy, H) + 1, mod(ox, W) + 1);
        [q, k] = max(sub(:));
        [ry, rx] = ind2sub(size(sub), k);
        dy = oy(ry); dx = ox(rx);
        % a maximum pinned to the edge of the window is not a peak, it is the
        % nearest point to a peak that lives outside - fall back to the prior
        onedge = (numel(ox) > 1 && (rx == 1 || rx == numel(ox))) || ...
                 (numel(oy) > 1 && (ry == 1 || ry == numel(oy)));
        if onedge
            dx = round(prior(1)); dy = round(prior(2)); ok = false;
        end
    end
    q = q / sqrt(sum(az(:).^2) * sum(bz(:).^2));   % normalized peak height
end

function w = hann_win(n)
    w = 0.5 * (1 - cos(2*pi*(0:n-1).'/(n-1)));
end

function s = rmfield_safe(s, f)
    if ischar(f), f = {f}; end
    for i = 1:numel(f)
        if isfield(s, f{i}), s = rmfield(s, f{i}); end
    end
end
