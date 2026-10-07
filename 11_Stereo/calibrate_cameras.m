function R = calibrate_cameras(opts)
% CALIBRATE_CAMERAS  Toolbox-free stereo-rig intrinsic calibration driver.
%
%   R = CALIBRATE_CAMERAS() detects the AMCalB dot target in every
%   calibration frame of both cameras, indexes each to the target grid,
%   runs Zhang's method per camera, prints the per-camera reprojection RMS
%   and intrinsics, and saves 9_Validation_DICe/cal_stereo.mat.
%
%   R = CALIBRATE_CAMERAS(opts) accepts:
%     .views    view numbers to use (default 0:65)
%     .cams     cameras to calibrate: 0, 1, or [0 1] (default [0 1])
%     .pitch    physical dot spacing (default 1 = unit pitch; intrinsics
%               are independent of it, so unit pitch is used until the real
%               AMCalB pitch is known)
%     .save     save cal_stereo.mat (default true)
%     .out      output .mat path (default 9_Validation_DICe/cal_stereo.mat)
%     .caldir   calibration image folder (default 9_Validation_DICe/calibration)
%     .stems    explicit list of file stems (cell): files are <stem>_<cam>.<ext>
%               in caldir; overrides .views and the generic listing, so a
%               driver can pick the calibration frames out of a mixed folder
%     .invert   invert images before detection (default false): white-dots-
%               on-black targets (Correlated Solutions / Vic-3D) become the
%               black-on-white pattern DETECT_DOTS expects
%     .ba       run the joint stereo bundle adjustment after the per-camera
%               solve (default true; needs both cameras). Set false to get
%               the raw per-camera Zhang result - useful for comparing the
%               two, not for measuring. See CALIBRATE_STEREO_BA.
%     .verbose  print progress (default true)
%
%   Output/saved struct fields: K0,K1 (intrinsics), dist0,dist1 ([k1 k2]),
%   cam(1).rvec/tvec and cam(2).rvec/tvec (per-view poses), kept0,kept1
%   (view numbers used per camera), rms0,rms1 (px), pitch, imsize. Each
%   cam(i) also carries detection statistics for quality gating:
%   n_views (candidate views on disk), n_indexed (survived detect+index),
%   n_failed_det, and iac_fallback (see CALIBRATE_ZHANG), plus views (the
%   per-view observations, aligned with rvec/tvec).
%
%   With two cameras the result is then refined jointly by
%   CALIBRATE_STEREO_BA and the refined intrinsics REPLACE the per-camera
%   ones; R.ba holds the full adjustment (including the relative pose, which
%   STEREO_GEOMETRY prefers over the per-view average) and R.rms_ba its
%   reprojection RMS. R.ba is [] if the adjustment was skipped or failed.
%
%   The 2D IC-GN engine is untouched; this is a standalone calibration
%   layer built on DETECT_DOTS + INDEX_GRID + CALIBRATE_ZHANG.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    if nargin < 1, opts = struct(); end
    vlist   = getf(opts, 'views',  0:65);
    pitch   = getf(opts, 'pitch',  1);
    dosave  = getf(opts, 'save',   true);
    caldir  = getf(opts, 'caldir', fullfile(root, '9_Validation_DICe', 'calibration'));
    outfile = getf(opts, 'out',    fullfile(root, '9_Validation_DICe', 'cal_stereo.mat'));
    verbose = getf(opts, 'verbose', true);
    cams    = getf(opts, 'cams', [0 1]);   % which cameras (0, 1, or both)
    stems   = getf(opts, 'stems', {});     % explicit stem list (mixed folders)
    invert  = getf(opts, 'invert', false); % white-on-black targets
    ba      = getf(opts, 'ba', true);      % joint stereo refinement (below)
    addpath(here);

    R = struct();
    imsize = [];
    for ci = 1:numel(cams)
        cam = cams(ci);
        V = struct('obj', {}, 'uv', {}, 'id', {});
        used = []; failed = [];
        files = {}; ids = [];
        if ~isempty(stems)                   % explicit stem list first
            exts = {'.tif', '.tiff', '.png', '.jpg', '.jpeg', '.bmp'};
            for q = 1:numel(stems)
                for e = 1:numel(exts)
                    f = fullfile(caldir, sprintf('%s_%d%s', stems{q}, cam, exts{e}));
                    if isfile(f)
                        files{end+1} = f; ids(end+1) = q; %#ok<AGROW>
                        break;
                    end
                end
            end
        end
        if isempty(files)
        for v = vlist                        % canonical AMCalB naming first
            f = fullfile(caldir, sprintf('AMCalB-%04d_%d.jpeg', v, cam));
            if isfile(f), files{end+1} = f; ids(end+1) = v; end %#ok<AGROW>
        end
        end
        if isempty(files)                    % any dot target: *_<cam>.<ext>
            dd = dir(caldir);
            names = sort({dd(~[dd.isdir]).name});
            pat = sprintf('_%d\\.(png|jpe?g|tiff?|bmp)$', cam);
            for q = 1:numel(names)
                if ~isempty(regexpi(names{q}, pat, 'once'))
                    files{end+1} = fullfile(caldir, names{q}); %#ok<AGROW>
                    ids(end+1) = numel(ids);                   %#ok<AGROW>
                end
            end
        end
        % Detection is embarrassingly parallel across views. parfor uses
        % all cores when the Parallel Computing Toolbox is present and
        % silently degrades to a plain serial loop when it is not, so the
        % toolbox-free contract is kept.
        nf = numel(files);
        objs = cell(1, nf); uvs = cell(1, nf); okv = false(1, nf);
        hs = zeros(1, nf); ws = zeros(1, nf); msgs = cell(1, nf);
        parfor vi = 1:nf
            try
                im = imread(files{vi});
            catch ME
                msgs{vi} = sprintf('cannot read %s (%s)', files{vi}, ME.message);
                continue;
            end
            hs(vi) = size(im, 1); ws(vi) = size(im, 2);
            if invert
                imd = double(im);
                im  = max(imd(:)) - imd;
            end
            D = detect_dots(im);
            Gv = index_grid(D, struct('pitch', pitch));
            if Gv.ok
                objs{vi} = Gv.obj; uvs{vi} = Gv.uv; okv(vi) = true;
            end
        end
        for vi = 1:nf
            v = ids(vi);
            if verbose && ~isempty(msgs{vi})
                fprintf('cam%d: %s - view skipped\n', cam, msgs{vi});
            end
            if okv(vi)
                V(end+1) = struct('obj', objs{vi}, 'uv', uvs{vi}, 'id', v); %#ok<AGROW>
                used(end+1) = v; %#ok<AGROW>
            else
                failed(end+1) = v; %#ok<AGROW>
            end
            if isempty(imsize) && hs(vi) > 0
                imsize = [hs(vi), ws(vi)];
            end
        end
        if verbose
            fprintf('cam%d: %d views indexed, %d failed detection %s\n', ...
                cam, numel(used), numel(failed), mat2str(failed));
        end

        C = calibrate_zhang(V, struct('imsize', imsize));
        kept_views = used(C.kept);
        drop_views = used(C.dropped);

        R.cam(ci).cam = cam;
        R.cam(ci).K = C.K; R.cam(ci).dist = C.dist;
        R.cam(ci).rvec = C.rvec; R.cam(ci).tvec = C.tvec;
        R.cam(ci).rms = C.rms; R.cam(ci).rms_view = C.rms_view;
        R.cam(ci).kept = kept_views; R.cam(ci).dropped = drop_views;
        % detection statistics for downstream quality gating: how many
        % candidate views existed on disk, how many survived detect+index,
        % and whether Zhang's closed form needed the focal fallback
        R.cam(ci).n_views      = nf;
        R.cam(ci).n_indexed    = numel(used);
        R.cam(ci).n_failed_det = numel(failed);
        R.cam(ci).iac_fallback = C.iac_fallback;
        % The observations behind the kept views, in the same order as
        % rvec/tvec. CALIBRATE_STEREO_BA needs them to refine both cameras
        % against one shared relative pose, and re-detecting the dots just
        % to get them back would cost as much as the calibration itself.
        R.cam(ci).views = V(C.kept);
        % Zhang's own per-view poses. The joint adjustment below replaces
        % rvec/tvec with poses that all share ONE relative pose, so they have
        % no spread left to show; these are independent per-view estimates,
        % and their spread is the rig-rigidity check of STEREO_EXTRINSICS.
        R.cam(ci).zhang = struct('rvec', {C.rvec}, 'tvec', {C.tvec}, 'kept', kept_views);

        if verbose
            fprintf(['cam%d: RMS %.4f px | fx %.1f fy %.1f | u0 %.1f v0 %.1f ' ...
                '| k1 %.4f k2 %.4f | %d/%d views | dropped-in-Zhang %s\n'], ...
                cam, C.rms, C.K(1,1), C.K(2,2), C.K(1,3), C.K(2,3), ...
                C.dist(1), C.dist(2), numel(kept_views), numel(used), mat2str(drop_views));
        end
    end

    R.K0 = R.cam(1).K; R.dist0 = R.cam(1).dist;
    R.kept0 = R.cam(1).kept; R.rms0 = R.cam(1).rms;
    if numel(cams) >= 2
        R.K1 = R.cam(2).K; R.dist1 = R.cam(2).dist;
        R.kept1 = R.cam(2).kept; R.rms1 = R.cam(2).rms;
    end
    R.pitch = pitch; R.imsize = imsize; R.cams = cams;

    % ---- joint stereo refinement ----
    % Per-camera Zhang is only the starting point: it cannot see the
    % focal/depth trade-off that a fronto-parallel target set leaves open.
    % Tying both cameras to one shared relative pose can. Failure here is
    % not fatal (per-camera result is returned) but it is reported.
    R.ba = [];
    if numel(cams) >= 2 && ba
        try
            R.ba = calibrate_stereo_ba(R, struct('verbose', verbose));
            R.K0 = R.ba.K0; R.dist0 = R.ba.dist0;
            R.K1 = R.ba.K1; R.dist1 = R.ba.dist1;
            R.cam(1).K = R.ba.K0; R.cam(1).dist = R.ba.dist0;
            R.cam(2).K = R.ba.K1; R.cam(2).dist = R.ba.dist1;
            R.rms_ba = R.ba.rms;
            % Per-view poses must match the adjusted K, otherwise rvec/tvec
            % users (STEREO_EXTRINSICS, re-runs on this output) get
            % inconsistent geometry.
            R = adopt_ba_poses(R);
        catch ME
            warning('calibrate_cameras:ba', ...
                'stereo bundle adjustment skipped (%s) - per-camera calibration kept.', ...
                ME.message);
        end
    end

    if dosave
        save(outfile, '-struct', 'R');
        if verbose, fprintf('saved %s\n', outfile); end
    end
end

function R = adopt_ba_poses(R)
% Replace the per-camera Zhang poses with the ones the joint adjustment
% actually solved, for the views it used. Camera 0 gets the view poses
% directly; camera 1 gets them composed with the shared relative pose.
% Views the adjustment did not see (not common to both cameras) are dropped
% rather than left stale.
    B = R.ba;
    for ci = 1:2
        keep = ismember([R.cam(ci).views.id], B.views);
        R.cam(ci).views = R.cam(ci).views(keep);
        R.cam(ci).kept  = R.cam(ci).kept(keep);
        if isfield(R.cam(ci), 'rms_view')
            R.cam(ci).rms_view = R.cam(ci).rms_view(keep);
        end
        n = numel(B.views);
        rv = cell(1, n); tv = cell(1, n);
        for k = 1:n
            Rw = rodr(B.rvec{k});  tw = B.tvec{k}(:);
            if ci == 2
                Rw = B.R * Rw;  tw = B.R * tw + B.T(:);
            end
            rv{k} = rotv(Rw);  tv{k} = tw;
        end
        R.cam(ci).rvec = rv;  R.cam(ci).tvec = tv;
    end
    R.kept0 = R.cam(1).kept;  R.kept1 = R.cam(2).kept;
    R.rms0  = R.ba.rms0;      R.rms1  = R.ba.rms1;
end

function M = rodr(v)
    v = double(v(:)); th = norm(v);
    if th < 1e-12, M = eye(3); return; end
    k = v/th; K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
    M = eye(3) + sin(th)*K + (1-cos(th))*K*K;
end

function v = rotv(M)
% Rotation vector of a rotation matrix, accurate at every angle. The angle
% comes from atan2(|w|, c) (acos loses digits near 0 and pi). Below 90 deg
% the axis is the skew part w = sin(th)*k; above it w shrinks to nothing as
% th -> pi (a board seen nearly upside down is th = pi - 1e-4) and dividing
% by sin(th) amplifies rounding, so the axis comes from the symmetric part
% (R + R')/2 - c*I = (1 - c)*k*k', with the sign taken from w.
    w = [M(3,2)-M(2,3); M(1,3)-M(3,1); M(2,1)-M(1,2)] / 2;
    c = (trace(M) - 1) / 2;
    th = atan2(norm(w), c);
    if th < 1e-12, v = [0;0;0]; return; end
    if c > 0
        v = th / norm(w) * w;
    else
        S = (M + M.') / 2 - c * eye(3);
        [~, i] = max(diag(S));
        k = S(:, i) / norm(S(:, i));
        if k.' * w < 0, k = -k; end
        v = th * k;
    end
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
