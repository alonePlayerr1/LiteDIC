function T = read_dice_gold(f)
% READ_DICE_GOLD  Parse a DICe solution/gold file into a column struct.
%
%   T = READ_DICE_GOLD(file) reads a DICe output file ('***' comment
%   lines, then a CSV whose first row names the columns) and returns a
%   struct with one field per column (names sanitized to valid MATLAB
%   identifiers), each an [n x 1] double vector. Columns are located by
%   NAME - never by position. Failed points are the caller's business
%   (DICe marks them MATCH = -1 / GAMMA = -1 in stereo output and
%   SIGMA = -1 in 2D output).
%
%   Used by VALIDATE_DICE_STEREO and the RUN_TESTS stereo-gold check.
%
%   See also VALIDATE_DICE_STEREO, STEREO_DIC.

    fid = fopen(f, 'r');
    if fid < 0, error('read_dice_gold:open', 'cannot open %s', f); end
    hdr = '';
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            fclose(fid);
            error('read_dice_gold:hdr', 'no header line in %s', f);
        end
        if ~isempty(ln) && ~strncmp(strtrim(ln), '***', 3)
            hdr = strtrim(ln); break;
        end
    end
    names = strsplit(hdr, ',');
    C = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', ',');
    fclose(fid);
    T = struct();
    for i = 1:numel(names)
        T.(matlab.lang.makeValidName(strtrim(names{i}))) = C{i};
    end
end
