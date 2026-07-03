function S = dic_scale(R, mm_per_px)
% DIC_SCALE  Express a 2D DIC result in millimetres instead of pixels.
%
%   S = DIC_SCALE(R, mm_per_px) takes a result struct from DIC_ENGINE or
%   DIC_ENGINE_BATCH and returns a copy whose length-dimension fields are
%   converted to millimetres by multiplying with mm_per_px:
%
%       displacements   S.U, S.V      (px -> mm)
%       grid coordinates S.gx, S.gy   (px -> mm)
%
%   Everything else is left exactly as it was: the correlation quality
%   (R.Q), the IC-GN iteration counts (R.NIT), and any strain fields
%   (exx, eyy, exy) that a caller has attached to the struct. Strains are
%   displacement gradients and therefore dimensionless - scaling the
%   displacements and the coordinates by the same factor leaves them
%   unchanged - so they pass through verbatim. Note the order: compute the
%   strains from the PIXEL result (PLS_STRAIN on R.U, R.V) BEFORE scaling
%   and attach them to R; this layer only carries pre-computed strain
%   fields through, it never differentiates. The scale factor and the unit
%   are recorded in S.mm_per_px and S.units for export/provenance.
%
%   mm_per_px is supplied by the user, typically from SCALE_FROM_LENGTH
%   (a known real-world distance between two picked points) or, later,
%   from a calibrated stereo baseline. It must be a positive finite scalar.
%
%   Back-compatibility: with mm_per_px = 1 every numeric field is returned
%   bit-for-bit identical to the input (multiplication by 1.0 is exact in
%   IEEE arithmetic), so adding this layer never changes the engine's
%   default pixel output. The IC-GN kernel is not involved at all - scaling
%   is a pure post-processing step on top of the untouched engine.
%
%   Example:
%       R = dic_engine(ref, cur, P);                        % px
%       [R.exx, R.eyy, R.exy] = pls_strain(R.U, R.V, ...);  % strains first
%       mm = scale_from_length(25.4, p1, p2);               % mm per px
%       S  = dic_scale(R, mm);   % U,V,gx,gy -> mm; strains carried through
%
%   See also DIC_ENGINE, DIC_ENGINE_BATCH, SCALE_FROM_LENGTH, PLS_STRAIN.

    if ~isstruct(R) || ~isscalar(R)
        error('dic_scale:input', 'R must be a scalar DIC result struct.');
    end
    if ~(isnumeric(mm_per_px) && isscalar(mm_per_px) && ...
            isfinite(mm_per_px) && mm_per_px > 0)
        error('dic_scale:scale', 'mm_per_px must be a positive finite scalar.');
    end

    S = R;                                   % copy; strains/Q/NIT pass through
    length_fields = {'U', 'V', 'gx', 'gy'};  % the only fields with px units
    for k = 1:numel(length_fields)
        f = length_fields{k};
        if isfield(S, f)
            S.(f) = S.(f) * mm_per_px;
        end
    end
    S.mm_per_px = mm_per_px;
    S.units = 'mm';
end
