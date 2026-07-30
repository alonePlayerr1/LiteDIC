function S = stereo_dic(files0, files1, G, P)
% STEREO_DIC  Stereo (3D) DIC pipeline layered on the untouched 2D engine.
%
%   S = STEREO_DIC(files0, files1, G, P) processes a synchronized
%   two-camera image series. files0/files1 are cell arrays of image paths
%   (camera 0 / camera 1), element 1 being the reference state. G is the
%   rig geometry from STEREO_GEOMETRY; P holds the DIC parameters (see
%   below). The 2D IC-GN engine (DIC_ENGINE) is used as-is for every
%   correlation; stereo adds only pre-alignment, composition, and
%   triangulation on top.
%
%   Per frame k the pipeline computes:
%     1. cross-camera disparity field D_k : cam0 (x,y) -> cam1
%        (STEREO_PAIR_MATCH on the synchronized pair k),
%     2. temporal displacement field u_k of camera 0 (DIC_ENGINE,
%        frame k vs frame 1), for k >= 2,
%     3. material-point triangulation on the reference grid:
%          reference:  X_ref = tri( x0,             x0 + D_1(x0) )
%          deformed :  X_k   = tri( x0 + u_k(x0),   x0 + u_k(x0)
%                                        + D_k(x0 + u_k(x0)) )
%        (D_k is evaluated AT the deformed position - the material point
%        moved), and the 3D displacement d3D = (X_k - X_ref) in mm via
%        G.pitch_mm.
%
%   Parameters (struct P) - passed to the engine:
%     P.subset_size, P.grid_step, P.search_range, P.black_thr,
%     P.zncc_min, P.max_iter, P.conv_tol
%   stereo-specific (optional):
%     P.search_temporal  integer search range for the temporal runs
%                        (default P.search_range)
%     P.zncc_min_cross   quality gate for cross-camera matching (default
%                        P.zncc_min; perspective between cameras lowers
%                        attainable ZNCC, so this is typically looser)
%     P.max_global       sanity clip for the global stereo shift, px
%
%   Output struct S:
%     S.gx, S.gy        reference grid (camera 0, from the cross match)
%     S.frames(k).D     cross-camera field (see STEREO_PAIR_MATCH)
%     S.frames(k).U0/V0 temporal cam0 displacement on S.gx/S.gy (k>=2)
%     S.frames(k).X_mm  triangulated 3D points, mm, [nY x nX x 3]
%     S.frames(k).d3D_mm 3D displacement vs reference, mm (k>=2)
%     S.frames(k).epi_rms  epipolar RMS of the matched pair, px
%     S.frames(k).epi_med  median epipolar residual (robust to clutter), px
%     S.frames(k).nvalid   valid matched nodes
%     S.X_ref_mm        reference-surface triangulation, mm
%     S.G, S.P          inputs, for provenance
%
%   See also STEREO_PAIR_MATCH, STEREO_GEOMETRY, TRIANGULATE_DLT, DIC_ENGINE.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(here, fullfile(root, '6_DIC'));

    nf = numel(files0);
    if numel(files1) ~= nf || nf < 1
        error('stereo_dic:files', 'files0/files1 must be equal-length, >= 1.');
    end

    Pt = P;                                  % temporal engine params
    if isfield(P, 'search_temporal') && ~isempty(P.search_temporal)
        Pt.search_range = P.search_temporal;
    end
    Pc = P;                                  % cross-camera engine params
    if isfield(P, 'zncc_min_cross') && ~isempty(P.zncc_min_cross)
        Pc.zncc_min = P.zncc_min_cross;
    end
    Pt = rmf(Pt, {'search_temporal', 'zncc_min_cross', 'max_global', ...
                  'shift_prior', 'shift_win'});
    Pc = rmf(Pc, {'search_temporal', 'zncc_min_cross'});

    pitch_mm = G.pitch_mm;
    if ~isfinite(pitch_mm)
        pitch_mm = 1;                        % fall back to pitch units
        warning('stereo_dic:pitch', 'pitch_mm unknown - output in pitch units.');
    end

    % ---- reference state ----
    I0r = load_gray(files0{1});
    I1r = load_gray(files1{1});
    S.frames = struct('D', {}, 'U0', {}, 'V0', {}, 'X_mm', {}, ...
        'd3D_mm', {}, 'epi_rms', {}, 'nvalid', {});

    Dr = stereo_pair_match(I0r, I1r, Pc);
    S.gx = Dr.gx; S.gy = Dr.gy;
    [Xg, Yg] = meshgrid(S.gx, S.gy);         % cam0 node coordinates

    [Xref, epi_r, epi_m, nval] = tri_field(Xg, Yg, Dr.U, Dr.V, G, pitch_mm);
    S.X_ref_mm = Xref;
    S.frames(1).D = Dr; S.frames(1).X_mm = Xref;
    S.frames(1).epi_rms = epi_r; S.frames(1).epi_med = epi_m;
    S.frames(1).nvalid = nval;
    S.frames(1).U0 = zeros(size(Xg)); S.frames(1).V0 = zeros(size(Xg));
    S.frames(1).d3D_mm = zeros(size(Xref));

    % ---- deformed states ----
    % Each frame's cross-camera pre-alignment is sought near the LAST
    % ACCEPTED one, starting from the reference pair's. The rig is rigid, so
    % that offset drifts (<= 25 px in dx, <= 1 px in dy between consecutive
    % frames, measured over 230 frames) but never jumps; see
    % STEREO_PAIR_MATCH for the measured failure mode this guards against.
    % Chaining tracks the drift, and because a frame with no credible peak
    % returns the prior unchanged (shift_ok false) instead of a boundary
    % artefact, a bad frame cannot poison the prior for the frames after it.
    % P.shift_win overrides the window.
    Pc.shift_prior = Dr.shift;
    nfallback = 0;
    for k = 2:nf
        I0 = load_gray(files0{k});
        I1 = load_gray(files1{k});

        Dk = stereo_pair_match(I0, I1, Pc);      % cross field of frame k
        if Dk.shift_ok
            Pc.shift_prior = Dk.shift;           % advance only on a real peak
        else
            nfallback = nfallback + 1;
        end
        Rk = dic_engine(I0r, I0, Pt);            % temporal, camera 0

        % temporal displacement interpolated onto the reference grid
        Fu = griddedInterpolant({Rk.gy, Rk.gx}, Rk.U, 'linear', 'none');
        Fv = griddedInterpolant({Rk.gy, Rk.gx}, Rk.V, 'linear', 'none');
        U0 = Fu(Yg, Xg); V0 = Fv(Yg, Xg);

        % disparity of frame k evaluated AT the deformed cam0 position
        Gu = griddedInterpolant({Dk.gy, Dk.gx}, Dk.U, 'linear', 'none');
        Gv = griddedInterpolant({Dk.gy, Dk.gx}, Dk.V, 'linear', 'none');
        p0x = Xg + U0;  p0y = Yg + V0;
        dU = Gu(p0y, p0x);  dV = Gv(p0y, p0x);

        [Xk, epi_r, epi_m, nval] = tri_field(p0x, p0y, dU, dV, G, pitch_mm);
        S.frames(k).D = Dk;
        S.frames(k).U0 = U0; S.frames(k).V0 = V0;
        S.frames(k).X_mm = Xk;
        S.frames(k).d3D_mm = Xk - Xref;
        S.frames(k).epi_rms = epi_r;
        S.frames(k).epi_med = epi_m;
        S.frames(k).nvalid = nval;
    end

    S.G = G; S.P = P;
    S.n_shift_fallback = nfallback;   % frames with no credible global peak
    if nfallback > 0
        fprintf(['stereo_dic: %d of %d frames had no credible global ' ...
            'pre-alignment peak inside the prior window - the previous ' ...
            'accepted shift was reused for those (see stereo_pair_match).\n'], ...
            nfallback, nf - 1);
    end
end

% ---------------------------------------------------------------------
function [Xmm, epi_rms, epi_med, nval] = tri_field(p0x, p0y, dU, dV, G, pitch_mm)
% Triangulate the correspondence (p0) -> (p0 + d) into [nY x nX x 3] mm.
    m = isfinite(p0x) & isfinite(p0y) & isfinite(dU) & isfinite(dV);
    uv0 = [p0x(m), p0y(m)];
    uv1 = [p0x(m) + dU(m), p0y(m) + dV(m)];
    nval = size(uv0, 1);
    Xmm = nan([size(p0x), 3]);
    if nval == 0, epi_rms = NaN; epi_med = NaN; return; end
    X = triangulate_dlt(uv0, uv1, G.K0, G.d0, G.K1, G.d1, G.R, G.T);
    for c = 1:3
        plane = nan(size(p0x));
        plane(m) = X(:, c) * pitch_mm;
        Xmm(:, :, c) = plane;
    end
    [epi_rms, epi_med] = epipolar_resid(uv0, uv1, G.F);
end

function [r, med] = epipolar_resid(uv0, uv1, F)
% RMS and median distance of uv1 to the epipolar lines of uv0 (px).
% The grid includes off-specimen clutter that passes the loose cross-camera
% quality gate, so the RMS is outlier-sensitive; the median is the robust
% field-quality number.
    n = size(uv0, 1);
    l = F * [uv0, ones(n,1)].';                      % lines in image 1
    num = sum([uv1, ones(n,1)].' .* l, 1);           % algebraic residual
    den = hypot(l(1,:), l(2,:));
    d = num ./ den;
    r = sqrt(mean(d.^2));
    med = median(abs(d));
end

function I = load_gray(f)
    I = imread(f);
    if isinteger(I), s = double(intmax(class(I))); else, s = 1; end
    I = double(I);
    if ndims(I) == 3, I = mean(I, 3); end
    I = I / s;
    if max(I(:)) > 1, I = I / max(I(:)); end   % non-integer inputs > 1
end

function s = rmf(s, ff)
    for i = 1:numel(ff)
        if isfield(s, ff{i}), s = rmfield(s, ff{i}); end
    end
end
