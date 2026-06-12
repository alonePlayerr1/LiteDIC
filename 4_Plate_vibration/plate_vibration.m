% Free vibration of a thin cantilever plate (Kirchhoff) - Rayleigh-Ritz.
% Boundary conditions: clamped at x = 0, all other edges free.
%
% The bending strain energy and the kinetic energy are integrated with
% Gauss-Legendre quadrature on an admissible polynomial basis (each basis
% function satisfies the clamped-edge conditions w = w' = 0 at x = 0).
%
% Sanity checks built into the output: for this narrow plate the first
% frequency must lie slightly ABOVE the Euler-Bernoulli cantilever beam
% value (a plate is stiffer than a beam due to Poisson coupling), and
% Rayleigh-Ritz converges to the true frequencies from above.
clear; close all; clc;

% Parameters
Lx = 1.2;         % length in x (clamped side), m
Ly = 0.2;         % length in y, m
h  = 0.01;        % plate thickness, m
E  = 7e10;        % Young's modulus, Pa
nu = 0.3;         % Poisson's ratio
rho = 2800;       % density, kg/m^3
D  = E * h^3 / (12 * (1 - nu^2));   % flexural rigidity

% Basis: X_i(x) = (x/Lx)^(i+1)  (clamped: w = w' = 0 at x = 0)
%        Y_j(y) = (y/Ly)^(j-1)  (free edges in y)
mx = 6;            % number of functions in x
ny = 4;            % number of functions in y
NB = mx * ny;

% Gauss-Legendre quadrature
ng = 30;
[xi, wq] = gauss_legendre(ng);
xn = 0.5 * Lx * (xi + 1); wxn = 0.5 * Lx * wq;   % nodes/weights in x
yn = 0.5 * Ly * (xi + 1); wyn = 0.5 * Ly * wq;   % nodes/weights in y
xn = xn(:)'; yn = yn(:)'; wxn = wxn(:)'; wyn = wyn(:)';

% Tabulated basis functions and derivatives at the quadrature nodes
Xv = zeros(mx, ng); X1 = zeros(mx, ng); X2 = zeros(mx, ng);
for i = 1:mx
    p = i + 1;
    Xv(i,:) = (xn/Lx).^p;
    X1(i,:) = p * xn.^(p-1) / Lx^p;
    X2(i,:) = p*(p-1) * xn.^(p-2) / Lx^p;
end
Yv = zeros(ny, ng); Y1 = zeros(ny, ng); Y2 = zeros(ny, ng);
for j = 1:ny
    q = j - 1;
    Yv(j,:) = (yn/Ly).^q;
    if q >= 1, Y1(j,:) = q * yn.^(q-1) / Ly^q; end
    if q >= 2, Y2(j,:) = q*(q-1) * yn.^(q-2) / Ly^q; end
end

% Fields of each basis function at the 2D nodes: w, w_xx, w_yy, w_xy
W   = zeros(ng, ng, NB);
WXX = zeros(ng, ng, NB);
WYY = zeros(ng, ng, NB);
WXY = zeros(ng, ng, NB);
for i = 1:mx
    for j = 1:ny
        k = (i-1)*ny + j;
        W(:,:,k)   = Xv(i,:).' * Yv(j,:);
        WXX(:,:,k) = X2(i,:).' * Yv(j,:);
        WYY(:,:,k) = Xv(i,:).' * Y2(j,:);
        WXY(:,:,k) = X1(i,:).' * Y1(j,:);
    end
end

% Stiffness K and mass M matrices (symmetric)
Wt = wxn.' * wyn;   % 2D quadrature weights (ng x ng)
K = zeros(NB, NB);
M = zeros(NB, NB);
for a = 1:NB
    for b = a:NB
        kern = WXX(:,:,a).*WXX(:,:,b) + WYY(:,:,a).*WYY(:,:,b) ...
             + nu*(WXX(:,:,a).*WYY(:,:,b) + WYY(:,:,a).*WXX(:,:,b)) ...
             + 2*(1-nu)*WXY(:,:,a).*WXY(:,:,b);
        K(a,b) = D * sum(sum(Wt .* kern));
        M(a,b) = rho * h * sum(sum(Wt .* (W(:,:,a).*W(:,:,b))));
        K(b,a) = K(a,b);
        M(b,a) = M(a,b);
    end
end

% Generalized eigenvalue problem
[V, Lam] = eig(K, M);
ev = real(diag(Lam));
[ev, si] = sort(ev, 'ascend');
V = V(:, si);
omega = sqrt(abs(ev));      % rad/s
f = omega / (2*pi);         % Hz

% First 5 natural frequencies
fprintf('Natural frequencies (Hz):\n');
for i = 1:min(5, NB)
    fprintf('f%d = %.3f Hz\n', i, f(i));
end

% Euler-Bernoulli cantilever beam comparison (narrow-plate check)
I_beam = Ly * h^3 / 12; A_beam = Ly * h;
f1_beam = (1.875^2/(2*pi)) * sqrt(E*I_beam/(rho*A_beam*Lx^4));
fprintf('Check: 1st cantilever beam frequency = %.3f Hz\n', f1_beam);

% First mode shape
[xg, yg] = meshgrid(linspace(0, Lx, 60), linspace(0, Ly, 30));
Wmode = zeros(size(xg));
for i = 1:mx
    for j = 1:ny
        k = (i-1)*ny + j;
        Wmode = Wmode + V(k,1) * (xg/Lx).^(i+1) .* (yg/Ly).^(j-1);
    end
end
Wmode = Wmode / max(abs(Wmode(:)));   % normalize

figure;
surf(xg, yg, Wmode);
xlabel('x (m)'); ylabel('y (m)'); zlabel('Normalized deflection w');
title(sprintf('First mode shape (f_1 = %.2f Hz)', f(1)));
shading interp; colorbar; view(-40, 30);

% ===== Local function (Gauss-Legendre nodes and weights on [-1,1]) =====
function [x, w] = gauss_legendre(n)
    k = 1:n-1;
    beta = k ./ sqrt(4*k.^2 - 1);
    T = diag(beta, 1) + diag(beta, -1);
    [Vv, Dd] = eig(T);
    [x, idx] = sort(diag(Dd));
    w = 2 * (Vv(1, idx).^2).';
end
