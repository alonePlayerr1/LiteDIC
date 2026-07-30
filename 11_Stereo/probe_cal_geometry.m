function probe_cal_geometry(dataDir, splitAt)
% PROBE_CAL_GEOMETRY  Diagnostic for a calibration set that does not fit
% a single K (uniform per-view residuals): per view and per camera, print
%   - lattice CHIRALITY: sign of cross(ex, ey) from INDEX_GRID. The donut
%     "L" is chiral, so all views of one physical target face must agree;
%     a sign flip = mirrored indexing (kills Zhang silently);
%   - fiducial arm ratio |long|/|short| (the L-geometry INDEX_GRID assumes);
%   - implied focal length from the view's own homography with the
%     principal point pinned at the image centre - if these scatter or
%     split into clusters, no single K can fit and the joint LM stalls at
%     a uniform residual floor;
%   - homography transfer error (indexing self-consistency).
% Also prints the per-view RMS stored in litedic_out\cal_stereo.mat, if
% present. Diagnostic only - modifies nothing.
%
%   probe_cal_geometry('C:\...\DIC 3D dataset', 400)

    if nargin < 2, splitAt = 400; end
    dd = dir(dataDir);
    names = {dd(~[dd.isdir]).name};
    stems = {};
    for q = 1:numel(names)
        t = regexpi(names{q}, '^(.*)_0\.(tif|tiff|png|jpe?g|bmp)$', 'tokens', 'once');
        if isempty(t) || strncmp(t{1}, '._', 2), continue; end
        n = regexp(t{1}, '(\d+)$', 'tokens', 'once');
        if ~isempty(n) && str2double(n{1}) < splitAt
            stems{end+1} = t{1}; %#ok<AGROW>
        end
    end
    fprintf('%d calibration stems\n', numel(stems));

    calfile = fullfile(dataDir, 'litedic_out', 'cal_stereo.mat');
    S = [];
    if isfile(calfile), S = load(calfile); end

    for cam = 0:1
        fprintf('\n===== camera %d =====\n', cam);
        fprintf('%-12s %5s %5s %8s %9s %9s %8s\n', 'stem', 'dots', 'chir', ...
            'arm-ratio', 'f_est1', 'f_est2', 'h-err');
        chir = []; fests = [];
        for q = 1:numel(stems)
            f = fullfile(dataDir, sprintf('%s_%d.tif', stems{q}, cam));
            if ~isfile(f), continue; end
            img = double(imread(f)); img = max(img(:)) - img;   % invert
            try
                D = detect_dots(img);
                G = index_grid(D);
            catch
                fprintf('%-12s  detect/index ERROR\n', stems{q});
                continue;
            end
            if ~G.ok
                fprintf('%-12s  grid not recovered\n', stems{q});
                continue;
            end
            c = sign(G.ex(1)*G.ey(2) - G.ex(2)*G.ey(1));
            % fiducial arm ratio in GRID indices (the physical L geometry)
            dk = find(D.is_donut);
            armr = NaN;
            if numel(dk) == 3
                dxy = [D.x(dk), D.y(dk)];
                dv = [norm(dxy(2,:)-dxy(1,:)), norm(dxy(3,:)-dxy(1,:)), ...
                      norm(dxy(3,:)-dxy(2,:))];
                dv = sort(dv);
                armr = dv(2) / dv(1);        % longer arm / shorter arm
            end
            H = dlt_h(G.obj(:,1:2), G.uv);
            herr = NaN; f1 = NaN; f2 = NaN;
            if ~isempty(H)
                pr = (H * [G.obj(:,1:2), ones(G.n,1)].').';
                pr = pr(:,1:2) ./ pr(:,3);
                herr = median(hypot(pr(:,1)-G.uv(:,1), pr(:,2)-G.uv(:,2)));
                [Himg, Wimg] = size(img);
                T = [1 0 -Wimg/2; 0 1 -Himg/2; 0 0 1];
                Hc = T * H;
                h1 = Hc(:,1); h2 = Hc(:,2);
                d = h1(3)*h2(3);
                if abs(d) > eps
                    v = -(h1(1)*h2(1) + h1(2)*h2(2))/d;
                    if v > 0, f1 = sqrt(v); end
                end
                d = h1(3)^2 - h2(3)^2;
                if abs(d) > eps
                    v = -((h1(1)^2+h1(2)^2) - (h2(1)^2+h2(2)^2))/d;
                    if v > 0, f2 = sqrt(v); end
                end
            end
            fprintf('%-12s %5d %5d %8.2f %9.0f %9.0f %8.2f\n', stems{q}, ...
                G.n, c, armr, f1, f2, herr);
            chir(end+1) = c; fests = [fests, f1, f2]; %#ok<AGROW>
        end
        fprintf(['cam%d summary: chirality +1: %d views, -1: %d views | ' ...
            'f median %.0f, IQR %.0f..%.0f (n=%d)\n'], cam, ...
            nnz(chir > 0), nnz(chir < 0), ...
            median(fests(isfinite(fests))), ...
            prctl(fests(isfinite(fests)), 25), prctl(fests(isfinite(fests)), 75), ...
            nnz(isfinite(fests)));
        if ~isempty(S)
            rv = S.cam(cam+1).rms_view;
            fprintf('cam%d stored per-view RMS: min %.2f | med %.2f | max %.2f px\n', ...
                cam, min(rv), median(rv), max(rv));
        end
    end
end

function H = dlt_h(XY, uv)
% Normalized-DLT homography (self-contained copy for this probe).
    n = size(XY, 1);
    if n < 4, H = []; return; end
    [X, Tx] = norm2(XY); [U, Tu] = norm2(uv);
    A = zeros(2*n, 9);
    for k = 1:n
        x = X(k,1); y = X(k,2); u = U(k,1); w = U(k,2);
        A(2*k-1, :) = [-x -y -1  0 0 0  u*x u*y u];
        A(2*k,   :) = [ 0 0 0  -x -y -1  w*x w*y w];
    end
    [~, ~, V] = svd(A, 0);
    H = Tu \ reshape(V(:, end), 3, 3).' * Tx;
    if abs(H(3,3)) > eps, H = H / H(3,3); end
end

function [Y, T] = norm2(X)
    c = mean(X, 1); d = X - c;
    s = sqrt(2) / mean(hypot(d(:,1), d(:,2)));
    T = [s 0 -s*c(1); 0 s -s*c(2); 0 0 1];
    Y = (T * [X, ones(size(X,1),1)].').';
    Y = Y(:, 1:2);
end

function v = prctl(x, p)
% percentile without toolboxes
    x = sort(x(:));
    if isempty(x), v = NaN; return; end
    v = x(max(1, min(numel(x), round(p/100 * numel(x)))));
end
