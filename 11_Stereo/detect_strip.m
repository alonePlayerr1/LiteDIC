function roi = detect_strip(img)
% DETECT_STRIP  Rectangular ROI around the textured specimen.
%
%   roi = DETECT_STRIP(img) takes a grayscale image (double) and returns
%   [x y w h]: a local-std texture map (decimated), then the strongest
%   contiguous run of textured columns, then the same over rows within those
%   columns. Works for a narrow speckled strip on a smooth background
%   regardless of orientation.
%
%   Split out of RUN_STEREO_DATASET (where it was a local function) so the
%   launcher can use the SAME specimen ROI when it measures the speckle to
%   size the subset. Measuring over the full frame instead let the estimate
%   land on background structure - see RESOLVE_RUN_PARAMS>SPECKLE_FRAME.
%
%   The threshold is referenced to the BACKGROUND level of the texture
%   profile, not to its mean. On a real specimen most of the profile IS
%   signal, so its mean sits well up the hump: measured on the ViC3D
%   sandwich panel the mean landed 60 % of the way from background to peak,
%   and speckle contrast tapers towards both ends of the panel, so a mean
%   threshold cut the ROI at y = 3020 of a 4112-row frame while Vic-3D
%   correlated down to y = 3937. A quarter of the specimen was not measured.
%
%   Each axis takes the highest-scoring run, not first..last above the
%   threshold (a bright frame edge at the top of these images pulled the
%   ROI to row 1). first..last was meant for a label or tape splitting the
%   texture profile mid-specimen; the low threshold covers that now: a dip
%   has to fall to within 10 % of background before it splits anything,
%   a label does not.
%
%   Base MATLAB, no toolboxes.
%
%   See also RUN_STEREO_DATASET, ESTIMATE_SPECKLE.
    [H, W] = size(img);
    dd = 4; w = 9;
    A = img(1:dd:end, 1:dd:end);
    K = ones(w) / w^2;
    mu = conv2(A, K, 'same');
    v  = max(conv2(A.^2, K, 'same') - mu.^2, 0);
    [c1, c2] = span_texture(mean(v, 1));
    [r1, r2] = span_texture(mean(v(:, c1:c2), 2).');
    pad = 3;
    x1 = max(1, (c1 - pad) * dd);  x2 = min(W, (c2 + pad) * dd);
    y1 = max(1, (r1 - pad) * dd);  y2 = min(H, (r2 + pad) * dd);
    roi = [x1, y1, x2 - x1, y2 - y1];
end

function [lo, hi] = span_texture(profile)
% Strongest contiguous run of a 1D texture profile above a threshold set
% 10 % of the way from the profile's background to its peak. Percentiles,
% not min/max: one saturated column at a frame edge would otherwise define
% "peak" and push the threshold up above the specimen itself.
    FRAC = 0.10;
    if isempty(profile), lo = 1; hi = 1; return; end
    p  = sort(profile);
    bg = p(max(1, round(0.05 * numel(p))));
    pk = p(max(1, round(0.98 * numel(p))));
    if pk <= bg, lo = 1; hi = numel(profile); return; end
    m = profile > bg + FRAC * (pk - bg);
    if ~any(m), lo = 1; hi = numel(profile); return; end
    d = diff([0, m, 0]);
    st = find(d == 1); en = find(d == -1) - 1;
    score = arrayfun(@(a, b) sum(profile(a:b)), st, en);
    [~, k] = max(score);
    lo = st(k); hi = en(k);
end
