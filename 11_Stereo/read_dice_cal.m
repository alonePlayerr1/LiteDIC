function C = read_dice_cal(xmlfile)
% READ_DICE_CAL  Parse a DICe legacy "vpml" stereo calibration file
% (tests/regression/.../cal.xml). Base MATLAB text parsing, no toolboxes.
%
%   C = READ_DICE_CAL(xmlfile) returns a struct with, per camera i in {1,2}:
%     C.cam(i).K       [3x3] intrinsics [fx 0 cx; 0 fy cy; 0 0 1]
%     C.cam(i).dist    [fs k1 k2 k3] as stored (DICe's own distortion slots)
%     C.cam(i).euler   [alpha beta gamma] orientation angles, degrees
%     C.cam(i).t       [tx ty tz] translation (mm)
%     C.cam(i).raw8    the 8 leading numbers as stored
%
%   The <camera> element stores eight numbers - cx cy fx fy fs k1 k2 k3 -
%   followed by an <orientation> of six - alpha beta gamma tx ty tz. This
%   reader only extracts them; the rotation convention is resolved and
%   checked separately (see STEREO_EXTRINSICS / the triangulation planarity
%   test), so nothing here assumes a particular Euler order.
%
%   See also STEREO_EXTRINSICS.

    txt = fileread(xmlfile);
    blocks = regexp(txt, '<camera\s+id="(\d+)">(.*?)</camera>', 'tokens');
    if numel(blocks) < 2
        error('read_dice_cal:parse', 'expected 2 <camera> blocks, found %d.', numel(blocks));
    end
    C = struct('file', xmlfile, 'cam', struct('K', {}, 'dist', {}, ...
        'euler', {}, 't', {}, 'raw8', {}));
    for b = 1:numel(blocks)
        inner = blocks{b}{2};
        ot = regexp(inner, '<orientation>(.*?)</orientation>', 'tokens', 'once');
        if isempty(ot)
            error('read_dice_cal:parse', 'no <orientation> in camera %d.', b);
        end
        pre = strtrim(regexprep(inner, '<orientation>.*', ''));   % text before <orientation>
        v8  = sscanf(pre, '%g');
        v6  = sscanf(ot{1}, '%g');
        if numel(v8) < 8 || numel(v6) < 6
            error('read_dice_cal:parse', 'camera %d: got %d + %d numbers.', b, numel(v8), numel(v6));
        end
        cx = v8(1); cy = v8(2); fx = v8(3); fy = v8(4);
        C.cam(b).K     = [fx 0 cx; 0 fy cy; 0 0 1];
        C.cam(b).dist  = v8(5:8).';
        C.cam(b).euler = v6(1:3).';
        C.cam(b).t     = v6(4:6).';
        C.cam(b).raw8  = v8(1:8).';
    end
end
