function [baseline_mm, R_rel, T_rel] = dice_stereo_baseline(C)
% DICE_STEREO_BASELINE  Camera0->camera1 relative pose and baseline (mm)
% from a parsed DICe stereo calibration.
%
%   [baseline_mm, R_rel, T_rel] = DICE_STEREO_BASELINE(C) takes the struct
%   from READ_DICE_CAL and returns the relative rotation R_rel, translation
%   T_rel (mm) and baseline |T_rel| (mm) between the two cameras.
%
%   Rotation convention: DICe's per-camera (alpha, beta, gamma) map to
%   R = Rz(gamma) * Ry(beta) * Rx(alpha), angles in degrees. This ordering
%   was pinned by matching the resulting relative rotation R_rel = R1*R0'
%   to the independently measured stereo rotation from STEREO_EXTRINSICS
%   (agreement < 0.002 in Frobenius norm; the other axis orders that share
%   it, being numerically identical here, give the same baseline). The
%   relative pose is world-frame-independent, so the recovered baseline is
%   unambiguous regardless of how DICe fixed its world origin.
%
%   See also READ_DICE_CAL, STEREO_EXTRINSICS.

    if isstruct(C) && isfield(C, 'file') && ~isfield(C, 'cam')
        C = read_dice_cal(C.file);
    elseif (ischar(C) || isstring(C))
        C = read_dice_cal(C);
    end
    R0 = eul_zyx(C.cam(1).euler); t0 = C.cam(1).t(:);
    R1 = eul_zyx(C.cam(2).euler); t1 = C.cam(2).t(:);
    R_rel = R1 * R0.';
    T_rel = t1 - R_rel * t0;
    baseline_mm = norm(T_rel);
end

function R = eul_zyx(a)
    R = rotz(a(3)) * roty(a(2)) * rotx(a(1));
end
function R = rotx(d), c = cosd(d); s = sind(d); R = [1 0 0; 0 c -s; 0 s c]; end
function R = roty(d), c = cosd(d); s = sind(d); R = [c 0 s; 0 1 0; -s 0 c]; end
function R = rotz(d), c = cosd(d); s = sind(d); R = [c -s 0; s c 0; 0 0 1]; end
