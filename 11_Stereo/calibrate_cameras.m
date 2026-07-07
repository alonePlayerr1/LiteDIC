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
%     .verbose  print progress (default true)
%
%   Output/saved struct fields: K0,K1 (intrinsics), dist0,dist1 ([k1 k2]),
%   cam(1).rvec/tvec and cam(2).rvec/tvec (per-view poses), kept0,kept1
%   (view numbers used per camera), rms0,rms1 (px), pitch, imsize.
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
    addpath(here);

    R = struct();
    imsize = [];
    for ci = 1:numel(cams)
        cam = cams(ci);
        V = struct('obj', {}, 'uv', {}, 'id', {});
        used = []; failed = [];
        for v = vlist
            f = fullfile(caldir, sprintf('AMCalB-%04d_%d.jpeg', v, cam));
            if ~isfile(f), continue; end
            im = imread(f);
            if isempty(imsize), imsize = [size(im,1), size(im,2)]; end
            D = detect_dots(im);
            G = index_grid(D, struct('pitch', pitch));
            if ~G.ok
                failed(end+1) = v; %#ok<AGROW>
                continue;
            end
            V(end+1) = struct('obj', G.obj, 'uv', G.uv, 'id', v); %#ok<AGROW>
            used(end+1) = v; %#ok<AGROW>
        end
        if verbose
            fprintf('cam%d: %d views indexed, %d failed detection %s\n', ...
                cam, numel(used), numel(failed), mat2str(failed));
        end

        C = calibrate_zhang(V);
        kept_views = used(C.kept);
        drop_views = used(C.dropped);

        R.cam(ci).cam = cam;
        R.cam(ci).K = C.K; R.cam(ci).dist = C.dist;
        R.cam(ci).rvec = C.rvec; R.cam(ci).tvec = C.tvec;
        R.cam(ci).rms = C.rms; R.cam(ci).rms_view = C.rms_view;
        R.cam(ci).kept = kept_views; R.cam(ci).dropped = drop_views;

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

    if dosave
        save(outfile, '-struct', 'R');
        if verbose, fprintf('saved %s\n', outfile); end
    end
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
