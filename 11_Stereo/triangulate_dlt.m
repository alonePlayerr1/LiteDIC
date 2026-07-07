function X = triangulate_dlt(uv0, uv1, K0, d0, K1, d1, R, T)
% TRIANGULATE_DLT  Linear (DLT) triangulation of matched stereo points into
% the camera-0 frame. Base MATLAB, no toolboxes.
%
%   X = TRIANGULATE_DLT(uv0, uv1, K0, d0, K1, d1, R, T) triangulates the
%   Nx2 matched image points uv0 (camera 0) and uv1 (camera 1). K0/K1 are
%   intrinsics, d0/d1 = [k1 k2] radial distortion, and (R, T) is the
%   camera-0 -> camera-1 relative pose. Points are undistorted and mapped
%   to normalized coordinates first, so the two projection matrices are the
%   metric P0 = [I|0] and P1 = [R|T]; X (Nx3) is in the same length units as
%   T (pitch units here).
%
%   See also STEREO_EXTRINSICS, READ_DICE_CAL.

    n0 = undistort_norm(uv0, K0, d0);
    n1 = undistort_norm(uv1, K1, d1);
    P0 = [eye(3), zeros(3,1)];
    P1 = [R, T(:)];
    N = size(uv0, 1);
    X = zeros(N, 3);
    for k = 1:N
        A = [ n0(k,1)*P0(3,:) - P0(1,:);
              n0(k,2)*P0(3,:) - P0(2,:);
              n1(k,1)*P1(3,:) - P1(1,:);
              n1(k,2)*P1(3,:) - P1(2,:) ];
        [~, ~, V] = svd(A, 0);
        h = V(:, end);
        X(k, :) = (h(1:3) / h(4)).';
    end
end

function xn = undistort_norm(uv, K, d)
% Pixel -> undistorted normalized coordinates (radial k1,k2 inverted by
% fixed-point iteration; distortion here is ~0.1%, converges in a few steps).
    x = (uv(:,1) - K(1,3)) / K(1,1);
    y = (uv(:,2) - K(2,3)) / K(2,2);
    k1 = d(1); k2 = d(2);
    xu = x; yu = y;
    for it = 1:6
        r2 = xu.^2 + yu.^2;
        f = 1 + k1*r2 + k2*r2.^2;
        xu = x ./ f; yu = y ./ f;
    end
    xn = [xu, yu];
end
