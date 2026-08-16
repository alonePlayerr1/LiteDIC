function [roi, bands] = detect_strip(img)
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
%   Each axis spans every SUBSTANTIAL run, not first..last above the
%   threshold (a bright frame edge pulls the ROI to row 1) and not only the
%   best run (a specimen with two speckled zones and a bare span keeps only
%   the first zone; DIC 3D dataset: zones at rows 4..1792 and 2968..3844 of
%   4112, Vic-3D AOI polygons 20..1724 and 3002..3800).
%
%   A run therefore counts as specimen when it scores at least KEEP of the
%   best run, and the ROI spans the first to the last that qualifies. On the
%   frame above the two speckle zones score 100 % and 42 %, everything else
%   (frame edges, testing machine, fixtures) 7.4 % or less.
%
%   Before that, dips NARROW relative to the bands they separate are closed
%   (marks on the specimen, not its end). Depth of the dip does not tell
%   them apart (a dark seam drops as low as a bare span), width does:
%
%     StereoSample5, seam mid-specimen   17 samples vs bands of 138 and 216
%     DIC 3D dataset, bare span          293        vs bands of 448 and 220
%     DIC 3D dataset, frame edge          34        vs a 19-sample edge run
%
%   Gaps up to a quarter of the shorter neighbouring band are closed, that
%   separates all three with a factor of two margin. Without it StereoSample5
%   was cut to the middle third by two surface marks.
%
%   The bare span between two zones is inside the ROI; nodes there fail to
%   correlate and are dropped. No mask: P.roi is one rectangle everywhere.
%
%   [roi, bands] = DETECT_STRIP(img) also returns the textured row bands the
%   ROI spans, in image pixels: an N x 2 array, one row per band. N > 1 means
%   the ROI bridges bare specimen, the caller should print that.
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
    [r1, r2, rbands] = span_texture(mean(v(:, c1:c2), 2).');
    pad = 3;
    x1 = max(1, (c1 - pad) * dd);  x2 = min(W, (c2 + pad) * dd);
    y1 = max(1, (r1 - pad) * dd);  y2 = min(H, (r2 + pad) * dd);
    roi = [x1, y1, x2 - x1, y2 - y1];
    bands = min(max(rbands * dd, 1), H);
end

function [lo, hi, bands] = span_texture(profile)
% Span of the substantial contiguous runs of a 1D texture profile above a
% threshold set 10 % of the way from the profile's background to its peak.
% Percentiles, not min/max: one saturated column at a frame edge would
% otherwise define "peak" and push the threshold up above the specimen
% itself.
%
% A run's score is the texture it carries (the sum, not the length), so a
% thin bright edge cannot qualify however bright it is, and a long faint
% band of background cannot either. See the header for the measured
% separation between the two populations.
    FRAC = 0.10;    % threshold, as a fraction of background-to-peak
    CLOSE = 0.25;   % close a gap this short relative to the shorter band beside it
    KEEP = 0.25;    % a run is specimen at this fraction of the best run's score
    if isempty(profile), lo = 1; hi = 1; bands = [1 1]; return; end
    p  = sort(profile);
    bg = p(max(1, round(0.05 * numel(p))));
    pk = p(max(1, round(0.98 * numel(p))));
    if pk <= bg
        lo = 1; hi = numel(profile); bands = [lo hi]; return;
    end
    m = profile > bg + FRAC * (pk - bg);
    if ~any(m)
        lo = 1; hi = numel(profile); bands = [lo hi]; return;
    end
    d = diff([0, m, 0]);
    st = find(d == 1); en = find(d == -1) - 1;
    [st, en] = close_narrow_gaps(st, en, CLOSE);
    score = arrayfun(@(a, b) sum(profile(a:b)), st, en);
    keep  = score >= KEEP * max(score);
    bands = [st(keep).', en(keep).'];
    lo = min(st(keep)); hi = max(en(keep));
end

function [st, en] = close_narrow_gaps(st, en, CLOSE)
% Merge neighbouring runs whose gap is short relative to the SHORTER of
% them. Relative, not absolute: the same seam is a few pixels on a close-up
% and a few hundred on a wide field, while the one measurement that stays
% comparable is how the gap compares with what it interrupts. Repeated to
% fixpoint, because merging lengthens a band and can qualify the next gap.
    merged = true;
    while merged && numel(st) > 1
        merged = false;
        for k = 1:numel(st)-1
            gap = st(k+1) - en(k) - 1;
            shorter = min(en(k) - st(k) + 1, en(k+1) - st(k+1) + 1);
            if gap <= CLOSE * shorter
                en(k) = en(k+1);
                st(k+1) = []; en(k+1) = [];
                merged = true;
                break;
            end
        end
    end
end
