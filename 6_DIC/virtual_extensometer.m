function E = virtual_extensometer(gx, gy, U, V, gauges, varargin)
% VIRTUAL_EXTENSOMETER  Gauge-length readings across arbitrary point pairs.
%
%   E = VIRTUAL_EXTENSOMETER(gx, gy, U, V, gauges) places one or more
%   virtual extensometers on a DIC displacement field and returns, per
%   gauge, the reference gauge length, the current length, the axial
%   (engineering) strain, and the RELATIVE displacement of the two ends
%   decomposed into an OPENING along the gauge axis and a SHEAR across it.
%
%   The point of a virtual extensometer is that it needs speckle only at
%   its two ENDS, not in between: straddle a feature that DIC cannot track
%   directly - a weld bead, a crack, a bond line, a bare band - with one
%   end on each bank, and the reading is the relative motion of the banks,
%   i.e. the opening/closing of the untrackable zone as a single number.
%   The field itself is NaN there. It is also the standard way to compare DIC against a physical clip
%   gauge / extensometer, which measures exactly this.
%
%   Inputs
%     gx, gy   grid coordinate vectors, px (as returned by dic_engine:
%              R.gx are columns/X, R.gy are rows/Y). For an accumulated
%              series field pass gx = X0(1,:), gy = Y0(:,1).
%     U, V     displacement fields, px (x right / y down), NaN where
%              invalid. Either a single state, size [numel(gy) numel(gx)],
%              or a stack of S states, size [numel(gy) numel(gx) S] - e.g.
%              one slice per loading step - which yields a reading curve.
%     gauges   K-by-4 array, each row [x1 y1 x2 y2]: the two gauge-end
%              coordinates in REFERENCE-frame px (frame-1 / undeformed).
%
%   Options (name/value)
%     'Scale'       mm-per-px; lengths (L0, L, opening, shear, du, dv) are
%                   returned in mm. Strain is dimensionless either way.
%                   Default 1 (outputs stay in px).
%     'SnapRadius'  px; if a gauge end lands where the field is NaN (just
%                   off a valid speckle node), snap it to the nearest node
%                   with a finite displacement within this radius. Default
%                   1.5*grid step (inferred from gx). 0 disables snapping;
%                   Inf snaps to the nearest valid node however far.
%     'Names'       K-by-1 cellstr of gauge labels (echoed in E.name).
%
%   Output struct E (all length quantities in the chosen unit; K gauges,
%   S states)
%     E.L0        K-by-1  reference gauge length
%     E.L         K-by-S  current gauge length per state
%     E.strain    K-by-S  axial engineering strain (L - L0)/L0
%     E.opening   K-by-S  relative end displacement PROJECTED on the
%                         reference gauge axis (extension; +ve = opening)
%     E.shear     K-by-S  relative end displacement across the axis
%                         (transverse sliding; signed, +ve = right of the
%                         axis when walking end1 -> end2)
%     E.du, E.dv  K-by-S  raw relative end displacement components
%     E.valid     K-by-S  true where BOTH ends had a finite reading
%     E.p1, E.p2  K-by-2  gauge-end coordinates actually used, px (after
%                         any snapping) - draw these to show the gauge
%     E.axis      K-by-2  unit reference gauge axis (end1 -> end2)
%     E.name      K-by-1  labels (or '' )
%     E.scale     the mm/px used
%
%   Orientation-agnostic: opening/shear are defined relative to each
%   gauge's own axis, so a vertical, horizontal or skew gauge all read
%   correctly (see test_flip_invariance.m for the field-level argument).
%   Base MATLAB only - no toolboxes.

    %% -------- options --------
    O.Scale      = 1;
    O.SnapRadius = [];
    O.Names      = {};
    for a = 1:2:numel(varargin)
        O.(varargin{a}) = varargin{a+1};
    end
    gx = gx(:).';  gy = gy(:).';
    if size(gauges, 2) ~= 4
        error('virtual_extensometer: gauges must be K-by-4 [x1 y1 x2 y2].');
    end
    if isempty(O.SnapRadius)
        if numel(gx) > 1, gstep = median(diff(gx)); else, gstep = 1; end
        O.SnapRadius = 1.5 * gstep;
    end
    K = size(gauges, 1);
    S = size(U, 3);
    if ~isequal(size(V, 1), size(U, 1)) || ~isequal(size(V, 2), size(U, 2)) ...
            || size(V, 3) ~= S
        error('virtual_extensometer: U and V must have the same size.');
    end
    if size(U, 1) ~= numel(gy) || size(U, 2) ~= numel(gx)
        error(['virtual_extensometer: U/V must be [numel(gy) numel(gx) S]; ' ...
               'got %dx%dx%d for gy=%d, gx=%d.'], size(U,1), size(U,2), S, ...
               numel(gy), numel(gx));
    end

    [GX, GY] = meshgrid(gx, gy);

    %% -------- per gauge --------
    E.L0     = nan(K, 1);
    E.L      = nan(K, S);
    E.strain = nan(K, S);
    E.opening = nan(K, S);
    E.shear  = nan(K, S);
    E.du     = nan(K, S);
    E.dv     = nan(K, S);
    E.valid  = false(K, S);
    E.p1     = nan(K, 2);
    E.p2     = nan(K, 2);
    E.axis   = nan(K, 2);
    E.name   = repmat({''}, K, 1);
    E.scale  = O.Scale;
    if ~isempty(O.Names), E.name(1:min(K, numel(O.Names))) = O.Names(:); end

    for g = 1:K
        p1 = gauges(g, 1:2);
        p2 = gauges(g, 3:4);
        d0 = p2 - p1;
        L0 = hypot(d0(1), d0(2));
        if L0 < eps
            warning('virtual_extensometer:zeroLength', ...
                'gauge %d has zero reference length - skipped.', g);
            continue;
        end
        ax = d0 / L0;                 % unit axis, end1 -> end2
        axp = [-ax(2), ax(1)];        % left-hand normal to the axis
        E.p1(g, :) = p1;  E.p2(g, :) = p2;
        E.L0(g) = L0 * O.Scale;
        E.axis(g, :) = ax;

        for s = 1:S
            [u1, v1] = sample(GX, GY, U(:,:,s), V(:,:,s), p1, O.SnapRadius);
            [u2, v2] = sample(GX, GY, U(:,:,s), V(:,:,s), p2, O.SnapRadius);
            if ~(isfinite(u1) && isfinite(v1) && isfinite(u2) && isfinite(v2))
                continue;             % leave NaN, valid stays false
            end
            q1 = p1 + [u1, v1];       % deformed ends
            q2 = p2 + [u2, v2];
            dcur = q2 - q1;
            L = hypot(dcur(1), dcur(2));
            drel = [u2 - u1, v2 - v1];         % relative end displacement
            E.L(g, s)      = L * O.Scale;
            E.strain(g, s) = (L - L0) / L0;    % dimensionless
            E.opening(g, s) = (drel(1)*ax(1)  + drel(2)*ax(2))  * O.Scale;
            E.shear(g, s)   = (drel(1)*axp(1) + drel(2)*axp(2)) * O.Scale;
            E.du(g, s) = drel(1) * O.Scale;
            E.dv(g, s) = drel(2) * O.Scale;
            E.valid(g, s) = true;
        end
    end
end

% ===================== local functions =====================

function [u, v] = sample(GX, GY, U, V, pt, radius)
% Displacement at a reference point pt = [x y] by bilinear interpolation,
% with a nearest-valid-node fallback: if the interpolated value is NaN
% (the point sits on/next to a masked hole) snap to the closest grid node
% with a finite displacement within `radius` px. This lets a gauge end be
% placed by eye near a speckle bank without landing exactly on a node.
    gx = GX(1, :);  gy = GY(:, 1).';
    u = interp2(gx, gy, U, pt(1), pt(2), 'linear', NaN);
    v = interp2(gx, gy, V, pt(1), pt(2), 'linear', NaN);
    if isfinite(u) && isfinite(v), return; end
    if radius <= 0, return; end
    fin = isfinite(U) & isfinite(V);
    if ~any(fin(:)), return; end
    d2 = (GX - pt(1)).^2 + (GY - pt(2)).^2;
    d2(~fin) = inf;
    [dmin, k] = min(d2(:));
    if sqrt(dmin) <= radius
        u = U(k);  v = V(k);
    end
end
