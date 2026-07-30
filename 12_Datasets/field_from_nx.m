function F = field_from_nx(csvfile, opts)
% FIELD_FROM_NX  Turn an FE (NX / Simcenter Nastran) surface solution into a
% GEN_DATASET displacement field, so the DIC engine can be validated against
% a finite-element ground truth on a real engineering case (a plate with a
% hole in tension) rather than only against analytic or third-party fields.
%
%   F = FIELD_FROM_NX(csvfile, opts)
%
%   csvfile : a CSV with a header line and columns  x,y,u,v  -- the node
%             coordinates and the two IN-PLANE displacement components of the
%             specimen FRONT SURFACE (the face a camera would see), in
%             millimetres, as taken from the solved NX model. (For the bundled
%             plate-with-hole model this is nx_plate_hole_field.csv, extracted
%             from the front face z = const of the CHEXA mesh.)
%
%   Returns a struct F with function handles
%       F.u(xpx, ypx) , F.v(xpx, ypx)
%   that give the displacement IN PIXELS at image pixel coordinates. Hand F to
%   GEN_DATASET through opts.field_fun and it renders an EXACT reference /
%   deformed speckle pair whose dense ground truth IS the FE field. The engine
%   then measures that pair and is compared back to the FE truth (see
%   NX_VALIDATE).
%
%   The mapping specimen(mm) -> image(px):
%     * a rectangular window  [xc +- Wmm/2] x [yc +- Hmm/2]  of the specimen,
%       centred on the hole, is imaged at S pixels per millimetre;
%     * the rigid-body MEAN of (u,v) over the window is removed -- a constant
%       offset changes no strain, and it keeps the deformed features inside the
%       frame;
%     * the field is multiplied by AMPLIFY, because a linear-static solve at a
%       modest load gives sub-0.1 px displacement; amplification turns it into
%       a few-pixel test WITHOUT changing the strain pattern (the ground truth
%       is scaled by exactly the same factor, so the ground truth stays exact).
%
%   Interpolation is by scatteredInterpolant (linear, with nearest-neighbour
%   extrapolation) -- base MATLAB, no toolbox -- so the handles are defined
%   everywhere the warp iteration and the strain differencing sample, including
%   the small disc of the hole (filled by extrapolation; a first-pass
%   simplification, see NOTE below).
%
%   opts (struct, all fields optional):
%     .scale      pixels per millimetre                     (default 20)
%     .win_mm     [Wmm Hmm] imaged window size, mm          (default [25 50])
%     .center_mm  [xc yc] window centre, mm                 (default [0 0])
%     .amplify    displacement multiplier                   (default 300)
%     .hole_r_mm  hole radius, mm; 0 = no hole              (default 2.5)
%
%   F also carries, for the caller:
%     F.u, F.v    the pixel-space displacement handles (above)
%     F.imsize    [H W] px  -- pass straight to gen_dataset opts.imsize
%     F.hole      struct .center_px [cx cy] (0-based), .r_px, .r_mm -- the
%                 hole disc in image pixels, so a caller (NX_VALIDATE) can
%                 blacken it and drop hole/rim nodes from the comparison
%     F.info      struct of all parameters + realised peak displacement (px)
%
%   NOTE (hole): the speckle is currently rendered over the whole window,
%   including the hole disc, whose displacement there is a smooth
%   extrapolation rather than true material. The stress-concentration lobes
%   AROUND the hole -- where the FE nodes are dense and the test matters -- are
%   faithful. Masking the hole out of both the image and the comparison is a
%   worthwhile refinement but is deliberately left out of this first pass.
%
%   See also GEN_DATASET, NX_VALIDATE, DIC_ENGINE.

    if nargin < 2, opts = struct(); end
    S       = getf(opts, 'scale',     20);
    win_mm  = getf(opts, 'win_mm',    [25 50]);
    ctr_mm  = getf(opts, 'center_mm', [0 0]);
    amp     = getf(opts, 'amplify',   300);
    hole_mm = getf(opts, 'hole_r_mm', 2.5);

    Wmm = win_mm(1);  Hmm = win_mm(2);
    xc  = ctr_mm(1);  yc  = ctr_mm(2);

    % ---- read the FE surface field (mm) ----
    T = readmatrix(csvfile);            % base MATLAB; skips the header row
    if size(T, 2) < 4
        error('field_from_nx:cols', ...
            'CSV %s must have columns x,y,u,v; got %d column(s).', ...
            csvfile, size(T, 2));
    end
    x = T(:,1);  y = T(:,2);  u = T(:,3);  v = T(:,4);

    % ---- restrict to the imaged window and remove the rigid mean ----
    inwin = abs(x - xc) <= Wmm/2 + 1e-9 & abs(y - yc) <= Hmm/2 + 1e-9;
    if nnz(inwin) < 10
        error('field_from_nx:empty', ...
            'Only %d FE nodes fall inside the %g x %g mm window at (%g,%g).', ...
            nnz(inwin), Wmm, Hmm, xc, yc);
    end
    umean = mean(u(inwin));
    vmean = mean(v(inwin));

    % ---- scattered interpolants in specimen mm (built on ALL nodes so the
    %      window edges extrapolate cleanly) ----
    Fu = scatteredInterpolant(x, y, u - umean, 'linear', 'nearest');
    Fv = scatteredInterpolant(x, y, v - vmean, 'linear', 'nearest');

    % ---- image geometry ----
    W = round(Wmm * S);  H = round(Hmm * S);   % [cols rows]
    x0 = xc - Wmm/2;                            % specimen x at pixel col 0
    y0 = yc - Hmm/2;                            % specimen y at pixel row 0

    % pixel (xpx,ypx) -> specimen mm ; displacement mm -> px is x S again,
    % so the whole displacement gets S * amp.
    px2mm_x = @(xpx) x0 + xpx / S;
    px2mm_y = @(ypx) y0 + ypx / S;
    k = S * amp;
    F.u = @(xpx, ypx) k * Fu(px2mm_x(xpx), px2mm_y(ypx));
    F.v = @(xpx, ypx) k * Fv(px2mm_x(xpx), px2mm_y(ypx));

    F.imsize = [H W];

    % ---- hole geometry in image pixels (0-based, matching gen_dataset's
    %      X=0..W-1 / Y=0..H-1 convention) ----
    F.hole = struct('center_px', [(xc - x0) * S, (yc - y0) * S], ...
                    'r_px', hole_mm * S, 'r_mm', hole_mm);

    % ---- realised peak displacement, for the log ----
    [XX, YY] = meshgrid(0:W-1, 0:H-1);
    up = F.u(XX, YY);  vp = F.v(XX, YY);
    F.info = struct('scale_px_per_mm', S, 'win_mm', win_mm, ...
        'center_mm', ctr_mm, 'amplify', amp, 'imsize', [H W], ...
        'n_nodes_in_window', nnz(inwin), ...
        'peak_u_px', max(abs(up(:))), 'peak_v_px', max(abs(vp(:))), ...
        'mean_removed_mm', [umean vmean], 'source', csvfile);

    fprintf(['field_from_nx: %d nodes in %g x %g mm window, %d x %d px @ ', ...
        '%g px/mm, amplify %g -> peak |u| %.2f px, |v| %.2f px\n'], ...
        nnz(inwin), Wmm, Hmm, W, H, S, amp, F.info.peak_u_px, F.info.peak_v_px);
end

function val = getf(s, name, dflt)
    if isfield(s, name) && ~isempty(s.(name)), val = s.(name); else, val = dflt; end
end
