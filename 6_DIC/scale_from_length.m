function mm_per_px = scale_from_length(known_mm, p1, p2)
% SCALE_FROM_LENGTH  Millimetres-per-pixel from a known real-world length.
%
%   mm_per_px = SCALE_FROM_LENGTH(known_mm, p1, p2) returns the physical
%   scale of an image from two points p1, p2 (each [x y] in pixels) whose
%   real separation in the scene is known_mm millimetres:
%
%       mm_per_px = known_mm / |p2 - p1|
%
%   mm_per_px = SCALE_FROM_LENGTH(known_mm) picks the two points
%   interactively with GINPUT on the current axes (base MATLAB graphics -
%   no toolbox). Show the reference frame first, then call this and click
%   the two ends of the known length.
%
%   Feed the result to DIC_SCALE to express displacements and coordinates
%   in millimetres. The two points should be as far apart as possible for
%   the best scale accuracy (the relative error in mm_per_px is the pixel
%   click error divided by the pixel span).
%
%   Example (a 25.4 mm gauge that spans 512 px in the image):
%       mm_per_px = scale_from_length(25.4, [100 300], [612 300]);  % 0.0496
%
%   See also DIC_SCALE, GINPUT.

    if ~(isnumeric(known_mm) && isscalar(known_mm) && ...
            isfinite(known_mm) && known_mm > 0)
        error('scale_from_length:known', ...
            'known_mm must be a positive finite scalar (millimetres).');
    end

    if nargin < 3
        title('click two points a known distance apart');
        pts = ginput(2);
        if size(pts, 1) < 2
            error('scale_from_length:aborted', 'two points were not picked.');
        end
        p1 = pts(1, :); p2 = pts(2, :);
    end

    p1 = p1(:).'; p2 = p2(:).';
    if numel(p1) ~= 2 || numel(p2) ~= 2
        error('scale_from_length:points', ...
            'p1 and p2 must each be [x y] pixel coordinates.');
    end

    d_px = hypot(p2(1) - p1(1), p2(2) - p1(2));
    if d_px < eps
        error('scale_from_length:degenerate', ...
            'the two points coincide (zero pixel distance).');
    end
    mm_per_px = known_mm / d_px;
end
