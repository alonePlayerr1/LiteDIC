function G = index_grid(D, params)
% INDEX_GRID  Assign integer grid coordinates to detected calibration dots
% using the three donut fiducials to anchor the target frame.
%
%   G = INDEX_GRID(D) takes the output of DETECT_DOTS and returns, for the
%   dots that belong to the regular lattice, their integer grid indices and
%   the matching planar object points (Z = 0) at unit pitch.
%
%   Output struct G:
%     G.ok      true if a consistent grid was recovered (needs 3 donuts)
%     G.i, G.j  integer grid indices (i along the long fiducial arm = X,
%               j along the short arm = Y); may be negative
%     G.uv      [n x 2] image points (px) matched to the indices
%     G.obj     [n x 3] object points [i j 0] * pitch, same row order as uv
%     G.n       number of indexed lattice dots
%     G.origin  the corner-donut image position [x y]
%     G.ex,G.ey one-step lattice vectors (px) along +i and +j
%
%   params (optional):
%     .pitch    physical dot spacing (default 1 = unit pitch; intrinsics do
%               not depend on it, so unit pitch is used until the real
%               AMCalB pitch is known)
%     .tol      neighbour snap tolerance as a fraction of the dot spacing
%               (default 0.40)
%
%   Frame: the three donuts form an L. The right-angle vertex is the
%   origin; the longer arm is the +X (i) axis, the shorter arm the +Y (j)
%   axis. Because the arms have different lengths this labelling is
%   unambiguous and identical in every view, so the object frame is glued
%   to the physical target across the whole calibration set (and across
%   both cameras of a synchronized pair).
%
%   Indexing grows outward from the origin dot: each placed dot predicts
%   its four lattice neighbours one LOCAL step ahead and snaps to the
%   nearest unclaimed detected dot within tol. Predicting only one local
%   step keeps the prediction accurate under perspective foreshortening,
%   and the tolerance gate rejects stray (non-lattice) detections and
%   tolerates missing dots.
%
%   See also DETECT_DOTS, CALIBRATE_ZHANG.

    G = struct('ok', false, 'i', [], 'j', [], 'uv', [], 'obj', [], ...
        'n', 0, 'origin', [], 'ex', [], 'ey', []);
    if nargin < 2, params = struct(); end
    pitch = getf(params, 'pitch', 1);
    tolf  = getf(params, 'tol', 0.40);

    P = [D.x(:), D.y(:)];
    N = size(P, 1);
    don = D.donut_xy;
    if size(don, 1) ~= 3 || N < 12 || ~isfinite(D.spacing), return; end

    % --- target frame from the three donuts ---
    [O, arm_long, arm_short] = donut_frame(don);
    nL = max(1, round(norm(arm_long  - O) / D.spacing));
    nS = max(1, round(norm(arm_short - O) / D.spacing));
    ex = (arm_long  - O) / nL;       % one lattice step along +i (X)
    ey = (arm_short - O) / nS;       % one lattice step along +j (Y)

    % --- start at the detected dot nearest the corner donut ---
    [~, s] = min(hypot(P(:,1) - O(1), P(:,2) - O(2)));
    tol = tolf * D.spacing;

    assigned = false(N, 1);
    I = zeros(N, 1); J = zeros(N, 1);
    EX = zeros(N, 2); EY = zeros(N, 2);
    assigned(s) = true; EX(s,:) = ex; EY(s,:) = ey;
    queue = s;
    dirs = [1 0; -1 0; 0 1; 0 -1];

    while ~isempty(queue)
        c = queue(1); queue(1) = [];
        pc = P(c,:); exc = EX(c,:); eyc = EY(c,:);
        for k = 1:4
            di = dirs(k,1); dj = dirs(k,2);
            pred = pc + di*exc + dj*eyc;
            dd = hypot(P(:,1) - pred(1), P(:,2) - pred(2));
            dd(assigned) = inf;
            [dm, m] = min(dd);
            if dm <= tol
                assigned(m) = true;
                I(m) = I(c) + di; J(m) = J(c) + dj;
                nex = exc; ney = eyc;
                if di ~= 0, nex = (P(m,:) - pc) / di; end
                if dj ~= 0, ney = (P(m,:) - pc) / dj; end
                EX(m,:) = nex; EY(m,:) = ney;
                queue(end+1) = m; %#ok<AGROW>
            end
        end
    end

    idx = find(assigned);
    if numel(idx) < 12, return; end
    G.ok = true;
    G.i = I(idx); G.j = J(idx);
    G.uv = P(idx, :);
    G.obj = [G.i * pitch, G.j * pitch, zeros(numel(idx), 1)];
    G.n = numel(idx);
    G.origin = O; G.ex = ex; G.ey = ey;
end

% ===================== helpers =====================

function [O, arm_long, arm_short] = donut_frame(don)
% Corner = the donut whose two arms are most orthogonal (the L vertex).
% The longer arm is returned first.
    best = inf; O = don(1,:); rest = don(2:3,:);
    for c = 1:3
        others = don(setdiff(1:3, c), :);
        v1 = others(1,:) - don(c,:);
        v2 = others(2,:) - don(c,:);
        cosa = abs(dot(v1, v2) / (norm(v1) * norm(v2) + eps));
        if cosa < best, best = cosa; O = don(c,:); rest = others; end
    end
    if norm(rest(1,:) - O) >= norm(rest(2,:) - O)
        arm_long = rest(1,:); arm_short = rest(2,:);
    else
        arm_long = rest(2,:); arm_short = rest(1,:);
    end
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
