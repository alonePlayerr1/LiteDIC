function [exx, eyy, exy] = pls_strain(U, V, step, win_px, use_green)
% PLS_STRAIN  Strains by pointwise least squares (Pan et al. 2009).
%
%   [exx, eyy, exy] = PLS_STRAIN(U, V, step, win_px, use_green) fits the
%   displacements inside a strain window with a plane (least squares) and
%   differentiates the fit ANALYTICALLY, so noise is averaged BEFORE
%   differentiation. U, V are displacement fields on a regular grid with
%   spacing `step` px (as returned by dic_engine); win_px is the strain
%   window in PIXELS; use_green = true switches to Green-Lagrange strains
%   (large rotations).
%
% Design choices that matter near the mask boundary:
%   1) the window is specified in PIXELS (win_px) and converted to grid
%      nodes -> refining the grid keeps the physical window size, so the
%      strain field is nearly independent of grid_step. Not
%      exactly: h is rounded to a whole number of nodes, and the fit is a
%      discrete sum over those nodes, so a coarser grid integrates the
%      same physical window with fewer samples. Measured in RUN_TESTS
%      (test_grid_refine): 2.7 % change in recovered strain amplitude
%      between a 4 px and a 10 px grid at a fixed 40 px window on a
%      120 px sinusoid, against 14.4 % if the window is fixed in NODES;
%   2) strains are computed ONLY at nodes with a valid displacement -
%      extrapolating the plane into the hole would fake large strain peaks
%      at the edge;
%   3) at the mask edge the fit requires enough points (>= win^2/3) and
%      coverage of the window along BOTH axes - guards against a
%      degenerate least-squares fit.
    [nY, nX] = size(U);
    h = max(1, round(win_px / (2*step)));        % half-window in nodes
    win = 2*h + 1;
    exx = nan(nY, nX); eyy = nan(nY, nX); exy = nan(nY, nX);
    valid = isfinite(U) & isfinite(V);
    need = max(6, ceil(win^2 / 3));
    [jf, iff] = meshgrid(-h:h, -h:h);            % jf = x (cols), iff = y (rows)
    Afull = pinv([ones(win^2,1), jf(:)*step, iff(:)*step]);  % partial-window path

    % ---- fast path: full interior windows via separable correlation ----
    % On the full symmetric window the plane-fit regressors (1, x, y) are
    % orthogonal, so the LSQ slopes reduce to slope_x = sum(x.*q)/sum(x.^2)
    % and slope_y = sum(y.*q)/sum(y.^2) - EXACTLY the Afull result, computed
    % for all points at once by two separable correlations. Applied only
    % where the whole window is valid and inside the array (the count test
    % below guarantees both), so NaNs never contaminate the fast results.
    cnt = conv2(ones(win,1), ones(1,win), double(valid), 'same');
    full_ok = valid & (cnt == win^2);
    if any(full_ok(:))
        xr  = (-h:h) * step;                     % window offsets, px
        Sq  = win * sum(xr.^2);                  % sum of x^2 (= sum of y^2)
        % correlation kernels are separable; conv2 flips them, so pass
        % flipped vectors (the odd factor changes sign, the even one not)
        kx_c = ones(win,1);        kx_r = flip(xr) / Sq;    % d/dx kernel
        ky_c = flip(xr).' / Sq;    ky_r = ones(1,win);      % d/dy kernel
        ux_f = conv2(kx_c, kx_r, U, 'same');
        uy_f = conv2(ky_c, ky_r, U, 'same');
        vx_f = conv2(kx_c, kx_r, V, 'same');
        vy_f = conv2(ky_c, ky_r, V, 'same');
        if use_green
            exx_f = ux_f + 0.5*(ux_f.^2 + vx_f.^2);
            eyy_f = vy_f + 0.5*(uy_f.^2 + vy_f.^2);
            exy_f = 0.5*(uy_f + vx_f + ux_f.*uy_f + vx_f.*vy_f);
        else
            exx_f = ux_f; eyy_f = vy_f; exy_f = 0.5*(uy_f + vx_f);
        end
        exx(full_ok) = exx_f(full_ok);
        eyy(full_ok) = eyy_f(full_ok);
        exy(full_ok) = exy_f(full_ok);
    end

    % ---- remaining points (mask edges, partial windows): original loop ----
    [ri, rj] = find(valid & ~full_ok);
    for k = 1:numel(ri)
        i = ri(k); j = rj(k);
        i0 = max(1, i-h); i1 = min(nY, i+h);
        j0 = max(1, j-h); j1 = min(nX, j+h);
        Ub = U(i0:i1, j0:j1); Vb = V(i0:i1, j0:j1);
        m = isfinite(Ub) & isfinite(Vb);
        if nnz(m) < need, continue; end
        if all(m(:)) && (i1-i0 == 2*h) && (j1-j0 == 2*h)
            au = Afull * Ub(:);              % full window: one multiply
            av = Afull * Vb(:);
        else                                 % mask edge: LSQ on valid points
            [jj, ii] = meshgrid((j0:j1)-j, (i0:i1)-i);
            xs = jj(m)*step; ys = ii(m)*step;
            if std(xs) < 0.5*step || std(ys) < 0.5*step, continue; end
            Xs = [ones(nnz(m),1), xs, ys];
            au = Xs \ Ub(m);
            av = Xs \ Vb(m);
        end
        ux = au(2); uy = au(3); vx = av(2); vy = av(3);
        if use_green
            F = [1+ux, uy; vx, 1+vy];
            E = 0.5 * (F' * F - eye(2));
            exx(i,j) = E(1,1); eyy(i,j) = E(2,2); exy(i,j) = E(1,2);
        else
            exx(i,j) = ux; eyy(i,j) = vy; exy(i,j) = 0.5*(uy + vx);
        end
    end
end
