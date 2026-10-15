function E = surf_strain(Xref, Xdef, opts)
% SURF_STRAIN  Surface (membrane) strain field from stereo DIC 3D points.
%
%   E = SURF_STRAIN(Xref, Xdef) computes the in-surface Green-Lagrange
%   strain field from two triangulated point grids [nY x nX x 3] (mm):
%   the reference surface Xref and a deformed state Xdef, both on the
%   same material grid (S.X_ref_mm and S.frames(k).X_mm of STEREO_DIC).
%
%   Method (pointwise least squares on the surface, the standard local
%   approach of commercial 3D DIC):
%     1. neighbours of each node in a (2*win+1)^2 grid window;
%     2. plane fit (SVD) to the REFERENCE neighbours -> local normal n;
%     3. tangent basis: e1 = global X axis projected onto the plane
%        (fallback: global Y when the plane is orthogonal to X),
%        e2 = n x e1 - so the axes are consistent across the surface;
%     4. 2D reference coordinates a_j of the neighbours in (e1,e2);
%     5. affine least squares  Xdef_j ~ t + F*a_j  with F [3x2] - the
%        surface deformation gradient;
%     6. C = F'*F and Green-Lagrange E = (C - I)/2 in the local frame:
%        E.exx (along e1), E.eyy (along e2), E.exy.
%   Rigid-body motion of an arbitrarily curved surface gives exactly
%   zero strain (F stays orthonormal); see the regression tests.
%
%   E = SURF_STRAIN(Xref, Xdef, opts) accepts:
%     .win      half-window in grid nodes (default 2 -> 5x5 neighbours;
%               larger = smoother, lower noise, less spatial resolution).
%               The default is in NODES, so its size on the specimen
%               follows the grid step: 5x5 nodes is 20 px at step 5 and
%               60 px at step 15. For a window of about 100 px whatever the
%               step (close to Vic-3D's strain filter) pass
%               win = round((100/gridStep - 1)/2), e.g. 4 (9x9) at step 11.
%     .min_pts  minimum finite neighbours for a solve (default 8)
%
%   Output struct E (all [nY x nX], NaN where not computable):
%     E.exx, E.eyy, E.exy   Green-Lagrange strains in the local frame
%     E.e1, E.e2            principal strains (e1 >= e2)
%     E.gmax                maximum shear, (e1 - e2)/2
%     E.win                 the window used
%
%   See also STEREO_DIC, PLS_STRAIN, RUN_STEREO_DATASET.

    if nargin < 3, opts = struct(); end
    win     = getf(opts, 'win', 2);
    min_pts = getf(opts, 'min_pts', 8);

    [nY, nX, three] = size(Xref);
    if three ~= 3 || ~isequal(size(Xdef), size(Xref))
        error('surf_strain:size', 'Xref/Xdef must be equal [nY x nX x 3] arrays.');
    end

    exx = nan(nY, nX); eyy = nan(nY, nX); exy = nan(nY, nX);
    finite = all(isfinite(Xref), 3) & all(isfinite(Xdef), 3);
    Pr = reshape(Xref, [], 3); Qd = reshape(Xdef, [], 3);

    for i = 1:nY
        for j = 1:nX
            if ~finite(i, j), continue; end
            ii = max(1, i-win):min(nY, i+win);
            jj = max(1, j-win):min(nX, j+win);
            [JJ, II] = meshgrid(jj, ii);
            idx = II(:) + (JJ(:) - 1) * nY;          % linear indices
            idx = idx(finite(idx));
            if numel(idx) < min_pts, continue; end

            P = Pr(idx, :);                          % reference neighbours
            Q = Qd(idx, :);                          % deformed neighbours
            Pc = mean(P, 1);
            Pd = P - Pc;

            [~, Sv, Vv] = svd(Pd, 'econ');
            sv = diag(Sv);
            if sv(2) < 1e-9 * max(sv(1), realmin)    % degenerate spread
                continue;
            end
            n = Vv(:, 3);                            % local plane normal
            % The sign of a singular vector is arbitrary (LAPACK picks it),
            % and e2 = n x e1 below flips with n, so exy would flip sign
            % from node to node. Fix: n_z > 0 (cameras look along +Z), or
            % the larger of |n_x|,
            % |n_y| for a surface seen edge-on. Then e1 ~ +X and e2 ~ +Y,
            % the same axes and shear sign as 2D DIC.
            if abs(n(3)) > 1e-6
                key = n(3);
            elseif abs(n(1)) >= abs(n(2))
                key = n(1);
            else
                key = n(2);
            end
            if key < 0, n = -n; end

            e1 = [1; 0; 0] - n(1) * n;               % global X onto plane
            if norm(e1) < 1e-6
                e1 = [0; 1; 0] - n(2) * n;
            end
            e1 = e1 / norm(e1);
            e2 = cross(n, e1);

            a = [Pd * e1, Pd * e2];                  % N x 2 tangent coords
            M = [a, ones(size(a, 1), 1)];            % N x 3 affine design
            C = M \ Q;                               % 3 x 3 LSQ solve
            F = C(1:2, :).';                         % 3 x 2 def. gradient
            C2 = F.' * F;
            exx(i, j) = 0.5 * (C2(1,1) - 1);
            eyy(i, j) = 0.5 * (C2(2,2) - 1);
            exy(i, j) = 0.5 *  C2(1,2);
        end
    end

    % principal strains and max shear (per node, closed form for 2x2 sym)
    mid = (exx + eyy) / 2;
    rad = sqrt(((exx - eyy) / 2).^2 + exy.^2);
    E = struct('exx', exx, 'eyy', eyy, 'exy', exy, ...
               'e1', mid + rad, 'e2', mid - rad, 'gmax', rad, 'win', win);
end

function v = getf(s, name, dflt)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = dflt; end
end
