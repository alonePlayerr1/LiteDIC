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
%      P.max_global (clip for the global shift sanity check, default 300 px).
%
%   Output struct D:
%     D.U, D.V    total disparity components, px (NaN where invalid)
%     D.gx, D.gy  camera-0 grid coordinates (same convention as DIC_ENGINE)
%     D.Q, D.NIT  engine quality / iteration fields
%     D.shift     the global integer [dx dy] that was pre-applied
%     D.peakq     normalized height of the correlation peak (diagnostic)
%
%   See also DIC_ENGINE, STEREO_DIC, STEREO_GEOMETRY.

    if nargin < 3, error('stereo_pair_match:args', 'img0, img1, P required.'); end
    maxg = 300;
    if isfield(P, 'max_global') && ~isempty(P.max_global), maxg = P.max_global; end
    P = rmfield_safe(P, 'max_global');   % dic_engine gets only its own fields

    % ---- 1) global integer offset by FFT cross-correlation ----
    [dx, dy, peakq] = global_shift(img0, img1);
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
    D.shift = [dx dy];
    D.peakq = peakq;
end

% ---------------------------------------------------------------------
function [dx, dy, q] = global_shift(a, b)
% Integer offset (dx, dy) with b(y+dy, x+dx) ~ a(y, x), from the peak of
% the circular FFT cross-correlation of the zero-mean, Hann-tapered
% images. c = ifft2(fft2(b) .* conj(fft2(a))) gives
%   c(1+dy, 1+dx) = sum_x b(x + [dx dy]) .* a(x)   (indices wrap),
% so the peak location IS the offset of b relative to a.
    [H, W] = size(a);
    wy = hann_win(H); wx = hann_win(W);
    Wn = wy * wx.';                       % separable taper, kills wrap edges
    az = (a - mean(a(:))) .* Wn;
    bz = (b - mean(b(:))) .* Wn;
    c = real(ifft2(fft2(bz) .* conj(fft2(az))));
    [q, k] = max(c(:));
    [iy, ix] = ind2sub([H, W], k);
    dy = iy - 1; dx = ix - 1;
    if dy > H/2, dy = dy - H; end         % unwrap to signed offsets
    if dx > W/2, dx = dx - W; end
    q = q / sqrt(sum(az(:).^2) * sum(bz(:).^2));   % normalized peak height
end

function w = hann_win(n)
    w = 0.5 * (1 - cos(2*pi*(0:n-1).'/(n-1)));
end

function s = rmfield_safe(s, f)
    if isfield(s, f), s = rmfield(s, f); end
end
