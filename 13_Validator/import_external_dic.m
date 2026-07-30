function E = import_external_dic(fname, opts)
% IMPORT_EXTERNAL_DIC  Read a displacement field produced by ANY DIC code.
%
%   E = IMPORT_EXTERNAL_DIC(fname) reads a text table (CSV / TSV /
%   semicolon / whitespace) with per-point results of an external DIC
%   algorithm (commercial package, open-source code or a neural network)
%   and returns a normalized struct for VALIDATE_EXTERNAL:
%     E.x, E.y   point coordinates, px
%     E.u, E.v   displacements, px (invalid points = NaN)
%     E.info     what was detected (delimiter, columns, n rows)
%
%   Column mapping:
%     * if the file has a header row, columns are found by NAME using a
%       tolerant alias list (x / x_px / coordinate_x / pos_x ...,
%       u / u_px / displacement_x / disp_x / dx ...), case-insensitive;
%     * without a header the first four numeric columns are taken as
%       x, y, u, v;
%     * or specify explicitly: opts.cols = [ix iy iu iv].
%
%   opts (all optional):
%     .cols       [ix iy iu iv] explicit column indices
%     .nan_below  u/v with |value| above this are invalid (default 1e6);
%                 also rows where a SIGMA-like column is negative are
%                 dropped automatically when such a column exists.
%
%   Comment lines starting with *, %, # are skipped (DICe writes ***).
%
%   See also VALIDATE_EXTERNAL, READMATRIX.

    if nargin < 2, opts = struct(); end
    cols_opt  = getf(opts, 'cols', []);
    nan_below = getf(opts, 'nan_below', 1e6);
    if ~isfile(fname), error('import_external_dic:file', 'File not found: %s', fname); end

    % ---- find the first non-comment line and sniff the delimiter ----
    fid = fopen(fname, 'r');
    n_skip = 0; first = '';
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        t = strtrim(line);
        if isempty(t) || startsWith(t, '*') || startsWith(t, '%') || startsWith(t, '#')
            n_skip = n_skip + 1; continue;
        end
        first = t; break;
    end
    fclose(fid);
    if isempty(first), error('import_external_dic:empty', 'No data rows in %s', fname); end

    if     contains(first, ',');  delim = ',';
    elseif contains(first, ';');  delim = ';';
    elseif contains(first, sprintf('\t')); delim = '\t';
    else,  delim = ' ';
    end

    % header row = the first line is not fully numeric
    toks = strtrim(strsplit(first, {',', ';', ' ', sprintf('\t')}));
    toks = toks(~cellfun(@isempty, toks));
    has_header = any(isnan(str2double(toks)) & ...
        ~ismember(lower(toks), {'nan', 'inf', '-inf'}));   % NaN literal = data

    hdr = {};
    if has_header
        hdr = strtrim(strsplit(first, delim));
        hdr = hdr(~cellfun(@isempty, hdr));
    end

    M = readmatrix(fname, 'FileType', 'text', 'NumHeaderLines', n_skip + has_header, ...
        'Delimiter', delim, 'ConsecutiveDelimitersRule', 'join');
    M = M(all(~isnan(M(:, 1:min(2, size(M, 2)))), 2), :);   % rows w/o coords are useless

    % ---- column mapping ----
    if ~isempty(cols_opt)
        ix = cols_opt(1); iy = cols_opt(2); iu = cols_opt(3); iv = cols_opt(4);
        how = 'explicit opts.cols';
    elseif has_header
        ix = find_alias(hdr, {'x', 'x_px', 'coordinate_x', 'pos_x', 'xc'});
        iy = find_alias(hdr, {'y', 'y_px', 'coordinate_y', 'pos_y', 'yc'});
        iu = find_alias(hdr, {'u', 'u_px', 'displacement_x', 'disp_x', 'dx', 'ux', 'u_mm'});
        iv = find_alias(hdr, {'v', 'v_px', 'displacement_y', 'disp_y', 'dy', 'uy', 'v_mm'});
        if any(cellfun(@isempty, {ix, iy, iu, iv}))
            error('import_external_dic:cols', ...
                ['Could not map columns by name in %s.\nHeader: %s\n' ...
                 'Pass opts.cols = [ix iy iu iv].'], fname, strjoin(hdr, ', '));
        end
        how = 'header names';
    else
        ix = 1; iy = 2; iu = 3; iv = 4;
        how = 'first four columns (no header)';
    end

    E.x = M(:, ix); E.y = M(:, iy); E.u = M(:, iu); E.v = M(:, iv);

    % ---- invalid points ----
    if has_header
        is = find_alias(hdr, {'sigma', 'quality', 'status'});
        if ~isempty(is)
            bad = M(:, is) < 0;
            E.u(bad) = NaN; E.v(bad) = NaN;
        end
    end
    big = abs(E.u) > nan_below | abs(E.v) > nan_below;
    E.u(big) = NaN; E.v(big) = NaN;

    E.info = struct('file', fname, 'delimiter', delim, 'mapping', how, ...
        'columns', [ix iy iu iv], 'n_points', numel(E.x), ...
        'n_valid', nnz(isfinite(E.u) & isfinite(E.v)));
end

function i = find_alias(hdr, names)
    i = [];
    h = lower(strtrim(hdr));
    for k = 1:numel(names)
        j = find(strcmpi(h, names{k}), 1);
        if ~isempty(j), i = j; return; end
    end
end

function v = getf(s, name, dflt)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = dflt; end
end
