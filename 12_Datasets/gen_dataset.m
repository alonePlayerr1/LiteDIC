function M = gen_dataset(opts)
% GEN_DATASET  Factory of synthetic DIC image pairs with EXACT ground truth.
%
%   M = GEN_DATASET(opts) renders speckle image pairs (reference +
%   deformed) together with mathematically exact dense displacement and
%   strain fields, for training and testing DIC algorithms (including
%   DL-DIC networks). The key property is that NO interpolation is used
%   anywhere in the imaging chain, so the ground truth carries no S-bias
%   and no resampling artifacts:
%
%     * the speckle pattern is an ANALYTIC function - a sum of Gaussian
%       spots I0(p) = bg + sum a_k exp(-|p-c_k|^2 / (2 s_k^2)) - that can
%       be evaluated exactly at any real coordinate;
%     * the deformed frame is I1(x) = I0(y(x)), where y solves the exact
%       material equation y + u(y) = x by fixed-point iteration (residual
%       < 1e-9 px, reported in M.qa);
%     * the ground truth is the analytic field u evaluated at integer
%       pixel coordinates of the reference frame - exactly the quantity
%       a DIC code measures; strains are derivatives of the analytic
%       field (step-1e-3 central differences, accurate to ~1e-9).
%
%   The camera model is applied AFTER the exact geometry: optional linear
%   illumination change on the deformed frame, additive Gaussian sensor
%   noise, clipping and 8/16-bit quantization.
%
%   opts (struct, all fields optional):
%     .outdir       output folder (default <this folder>/output)
%     .n_pairs      number of image pairs (8)
%     .imsize       [H W] pixels ([512 512])
%     .seed         RNG seed; same seed => bit-identical dataset (1)
%     .field        'mixed' (default) | 'rigid' | 'affine' | 'sine' | 'gauss'
%     .max_disp     displacement cap, px (4)
%     .max_strain   strain cap (0.02)
%     .dots_per_kpx speckle density per 1000 px (8)
%     .sigma_rng    spot radius range, px ([1.0 2.2])
%     .amp_rng      spot amplitude range ([0.4 0.95])
%     .bg           background level (0.12)
%     .noise_std    Gaussian sensor noise std (0.005; 0 = clean)
%     .illum_grad   max relative illumination tilt on frame 2 (0.03)
%     .bits         PNG bit depth, 8 or 16 (16)
%     .csv          also write the GT as CSV per pair (false)
%     .field_fun    external displacement field: a struct with handles .u(x,y)
%                   and .v(x,y) in pixels (and optional .info), overriding the
%                   built-in .field library. Set opts.imsize to match. Used to
%                   drive the factory from an FE solution (see FIELD_FROM_NX),
%                   so the DIC engine can be validated against a finite-element
%                   ground truth. Absent/empty = built-in analytic fields.
%
%   Output per pair i (in outdir):
%     pair_%04d_ref.png / _def.png   the images
%     pair_%04d_gt.mat               U, V, EXX, EYY, EXY [H x W], field info
%     pair_%04d_gt.csv               (if opts.csv) x,y,u,v,exx,eyy,exy
%   plus dataset_meta.json (all parameters + per-pair QA) and M:
%     M.files     cell of written file names
%     M.qa(i)     .inv_residual (px), .max_disp, .max_strain of each pair
%
%   See also DIC_ENGINE, BIAS_SCURVE, RUN_TESTS.

    here = fileparts(mfilename('fullpath'));
    if nargin < 1, opts = struct(); end
    o.outdir       = getf(opts, 'outdir', fullfile(here, 'output'));
    o.n_pairs      = getf(opts, 'n_pairs', 8);
    o.imsize       = getf(opts, 'imsize', [512 512]);
    o.seed         = getf(opts, 'seed', 1);
    o.field        = getf(opts, 'field', 'mixed');
    o.max_disp     = getf(opts, 'max_disp', 4);
    o.max_strain   = getf(opts, 'max_strain', 0.02);
    o.dots_per_kpx = getf(opts, 'dots_per_kpx', 8);
    o.sigma_rng    = getf(opts, 'sigma_rng', [1.0 2.2]);
    o.amp_rng      = getf(opts, 'amp_rng', [0.4 0.95]);
    o.bg           = getf(opts, 'bg', 0.12);
    o.noise_std    = getf(opts, 'noise_std', 0.005);
    o.illum_grad   = getf(opts, 'illum_grad', 0.03);
    o.bits         = getf(opts, 'bits', 16);
    o.csv          = getf(opts, 'csv', false);
    o.field_fun    = getf(opts, 'field_fun', []);   % external field (see below)

    if ~isfolder(o.outdir), mkdir(o.outdir); end
    H = o.imsize(1); W = o.imsize(2);
    [X, Y] = meshgrid(0:W-1, 0:H-1);
    types = {'rigid', 'affine', 'sine', 'gauss'};

    M.files = {}; M.qa = struct('inv_residual', {}, 'max_disp', {}, 'max_strain', {});
    t_all = tic;
    for i = 1:o.n_pairs
        rng(o.seed * 100000 + i, 'twister');

        % ---- speckle set (analytic pattern) ----
        nd = round(o.dots_per_kpx * H * W / 1000);
        margin = 3;
        cx = margin + (W - 2*margin) * rand(nd, 1);
        cy = margin + (H - 2*margin) * rand(nd, 1);
        sg = o.sigma_rng(1) + diff(o.sigma_rng) * rand(nd, 1);
        am = o.amp_rng(1)   + diff(o.amp_rng)   * rand(nd, 1);

        % ---- displacement field ----
        % By default the field is an analytic handle from the built-in
        % library (make_field). If the caller passes opts.field_fun (a struct
        % with .u and .v function handles of pixel coordinates, e.g. an FE
        % field from FIELD_FROM_NX), that field is used instead -- same handle
        % contract, so the exact-warp / exact-strain machinery below is
        % unchanged. Without field_fun the built-in library is used.
        if ~isempty(o.field_fun)
            ufun  = o.field_fun.u;
            vfun  = o.field_fun.v;
            ftype = 'external';
            if isfield(o.field_fun, 'info'), finfo = o.field_fun.info;
            else, finfo = struct(); end
        else
            if strcmpi(o.field, 'mixed'), ftype = types{mod(i-1, numel(types)) + 1};
            else, ftype = lower(o.field); end
            [ufun, vfun, finfo] = make_field(ftype, W, H, o.max_disp, o.max_strain);
        end

        % ---- exact inverse warp: y + u(y) = x ----
        Yx = X; Yy = Y;
        for it = 1:25
            Yx = X - ufun(Yx, Yy);
            Yy = Y - vfun(Yx, Yy);
        end
        inv_res = max(max(hypot(Yx + ufun(Yx, Yy) - X, Yy + vfun(Yx, Yy) - Y)));

        % ---- render both frames (windowed 5-sigma stamping, exact) ----
        umax = max([max(abs(ufun(X, Y)), [], 'all'), ...
                    max(abs(vfun(X, Y)), [], 'all')]);
        I0 = stamp(X,  Y,  cx, cy, sg, am, o.bg, 0);
        I1 = stamp(Yx, Yy, cx, cy, sg, am, o.bg, umax);

        % ---- ground truth at reference pixel coordinates ----
        U = ufun(X, Y); V = vfun(X, Y);
        h = 1e-3;
        EXX = (ufun(X+h, Y) - ufun(X-h, Y)) / (2*h);
        EYY = (vfun(X, Y+h) - vfun(X, Y-h)) / (2*h);
        EXY = ((ufun(X, Y+h) - ufun(X, Y-h)) + (vfun(X+h, Y) - vfun(X-h, Y))) / (4*h);

        % ---- camera model (geometry is already final) ----
        if o.illum_grad > 0
            gxr = o.illum_grad * (2*rand - 1); gyr = o.illum_grad * (2*rand - 1);
            I1 = I1 .* (1 + gxr*(X/W - 0.5) + gyr*(Y/H - 0.5));
        end
        if o.noise_std > 0
            I0 = I0 + o.noise_std * randn(H, W);
            I1 = I1 + o.noise_std * randn(H, W);
        end
        I0 = min(max(I0, 0), 1); I1 = min(max(I1, 0), 1);

        % ---- write ----
        f_ref = fullfile(o.outdir, sprintf('pair_%04d_ref.png', i));
        f_def = fullfile(o.outdir, sprintf('pair_%04d_def.png', i));
        f_gt  = fullfile(o.outdir, sprintf('pair_%04d_gt.mat',  i));
        if o.bits == 8
            imwrite(uint8(round(I0 * 255)),  f_ref);
            imwrite(uint8(round(I1 * 255)),  f_def);
        else
            imwrite(uint16(round(I0 * 65535)), f_ref);
            imwrite(uint16(round(I1 * 65535)), f_def);
        end
        field_type = ftype; field_info = finfo; %#ok<NASGU>
        save(f_gt, 'U', 'V', 'EXX', 'EYY', 'EXY', 'field_type', 'field_info');
        M.files = [M.files, {f_ref, f_def, f_gt}];
        if o.csv
            f_csv = fullfile(o.outdir, sprintf('pair_%04d_gt.csv', i));
            write_gt_csv(f_csv, X, Y, U, V, EXX, EYY, EXY);
            M.files{end+1} = f_csv;
        end
        M.qa(i).inv_residual = inv_res;
        M.qa(i).max_disp     = umax;
        M.qa(i).max_strain   = max(max(abs([EXX(:); EYY(:); EXY(:)])));
        fprintf(['pair %3d/%d: %-6s  max|u| %.2f px, max|strain| %.4f, ' ...
            'inverse residual %.1e px\n'], i, o.n_pairs, ftype, ...
            M.qa(i).max_disp, M.qa(i).max_strain, inv_res);
    end

    meta = o; meta.qa = M.qa; meta.generated = char(datetime('now'));
    if isfield(meta, 'field_fun') && ~isempty(meta.field_fun)
        % function handles cannot be JSON-encoded: keep the .info summary
        % (if any) and drop the handles.
        if isstruct(meta.field_fun) && isfield(meta.field_fun, 'info')
            meta.field_fun = meta.field_fun.info;
        else
            meta.field_fun = 'external field handles (not serialized)';
        end
    end
    fid = fopen(fullfile(o.outdir, 'dataset_meta.json'), 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true)); fclose(fid);
    fprintf('Dataset: %d pair(s) in %.1f s -> %s\n', o.n_pairs, toc(t_all), o.outdir);
end

% =====================================================================
function I = stamp(Px, Py, cx, cy, sg, am, bg, pad_extra)
% Exact evaluation of the Gaussian-spot pattern at coordinates (Px,Py),
% accumulated per spot over a 5-sigma window (+ warp padding). Truncation
% tail < 4e-6 - far below the 16-bit quantization step.
    [H, W] = size(Px);
    I = bg * ones(H, W);
    for k = 1:numel(cx)
        r  = ceil(5 * sg(k) + pad_extra) + 2;
        x0 = max(1, round(cx(k)) - r); x1 = min(W, round(cx(k)) + r);
        y0 = max(1, round(cy(k)) - r); y1 = min(H, round(cy(k)) + r);
        px = Px(y0:y1, x0:x1) - cx(k);
        py = Py(y0:y1, x0:x1) - cy(k);
        I(y0:y1, x0:x1) = I(y0:y1, x0:x1) + ...
            am(k) * exp(-(px.^2 + py.^2) / (2 * sg(k)^2));
    end
end

function [ufun, vfun, info] = make_field(ftype, W, H, max_disp, max_strain)
% Analytic displacement-field library. All randomness comes from the
% caller's RNG state, so datasets are reproducible by seed.
    cxc = (W-1)/2; cyc = (H-1)/2;
    switch ftype
        case 'rigid'
            tx = max_disp * (2*rand - 1) * 0.7;
            ty = max_disp * (2*rand - 1) * 0.7;
            th = min(0.3 * max_disp / hypot(W, H) * 2, 2e-3) * (2*rand - 1);
            ufun = @(x, y) tx + (cos(th)-1).*(x-cxc) - sin(th).*(y-cyc);
            vfun = @(x, y) ty + sin(th).*(x-cxc) + (cos(th)-1).*(y-cyc);
            info = struct('tx', tx, 'ty', ty, 'theta', th);
        case 'affine'
            A = max_strain * (2*rand(2, 2) - 1);
            tx = max_disp * (2*rand - 1) * 0.5;
            ty = max_disp * (2*rand - 1) * 0.5;
            ufun = @(x, y) tx + A(1,1)*(x-cxc) + A(1,2)*(y-cyc);
            vfun = @(x, y) ty + A(2,1)*(x-cxc) + A(2,2)*(y-cyc);
            info = struct('A', A, 'tx', tx, 'ty', ty);
        case 'sine'
            lx = 24 + rand * 72;  ly = 24 + rand * 72;
            ax = min(max_disp, max_strain * lx / (2*pi)) * (0.5 + 0.5*rand);
            ay = min(max_disp, max_strain * ly / (2*pi)) * (0.5 + 0.5*rand);
            phx = 2*pi*rand; phy = 2*pi*rand;
            ufun = @(x, y) ax * sin(2*pi*x/lx + phx) .* cos(2*pi*y/(ly*1.3) + phy);
            vfun = @(x, y) ay * cos(2*pi*x/(lx*1.3) + phx) .* sin(2*pi*y/ly + phy);
            info = struct('lambda', [lx ly], 'amp', [ax ay], 'phase', [phx phy]);
        case 'gauss'
            nb = 2 + randi(2);
            bx = W * rand(nb, 1); by = H * rand(nb, 1);
            wd = 15 + 35 * rand(nb, 1);
            aa = min(max_disp, max_strain * wd / 0.607) .* (2*rand(nb, 1) - 1);
            dirs = 2*pi*rand(nb, 1);
            ufun = @(x, y) gsum(x, y, bx, by, wd, aa .* cos(dirs));
            vfun = @(x, y) gsum(x, y, bx, by, wd, aa .* sin(dirs));
            info = struct('centers', [bx by], 'width', wd, 'amp', aa, 'dir', dirs);
        otherwise
            error('gen_dataset:field', 'Unknown field type "%s".', ftype);
    end
end

function s = gsum(x, y, bx, by, wd, aa)
    s = zeros(size(x));
    for b = 1:numel(bx)
        s = s + aa(b) * exp(-((x-bx(b)).^2 + (y-by(b)).^2) / (2 * wd(b)^2));
    end
end

function write_gt_csv(fname, X, Y, U, V, EXX, EYY, EXY)
    fid = fopen(fname, 'w');
    if fid < 0, warning('gen_dataset:csv', 'Cannot write %s', fname); return; end
    fprintf(fid, 'x_px,y_px,u_px,v_px,exx,eyy,exy\n');
    T = [X(:), Y(:), U(:), V(:), EXX(:), EYY(:), EXY(:)];
    fprintf(fid, '%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g\n', T.');
    fclose(fid);
end

function v = getf(s, name, dflt)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = dflt; end
end
