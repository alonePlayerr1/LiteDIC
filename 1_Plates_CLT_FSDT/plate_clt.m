% Bending of a laminated composite plate - Classical Lamination Theory (CLT).
%
% Simply supported rectangular [0/90/0] laminate under uniform pressure,
% solved with a Navier double sine series. Outputs: ABD stiffness matrices,
% center deflection, curvatures, bending moments and ply-by-ply stresses
% in material axes.
clear all;
clc;

% Plate parameters
a = 1;            % plate length, m
b = 0.5;          % plate width, m
h = 0.006;        % total laminate thickness, m
h_k = 0.002;      % ply thickness, m
theta = [0, 90, 0];          % ply orientation angles, deg
n_layers = length(theta);    % number of plies
q0 = 10e3;        % uniform pressure, N/m^2

% Material properties (unidirectional ply)
E1 = 140e9;       % modulus along fibers, Pa
E2 = 10e9;        % modulus across fibers, Pa
G12 = 5e9;        % in-plane shear modulus, Pa
nu12 = 0.3;       % major Poisson's ratio
nu21 = nu12 * E2 / E1;

% Ply stiffness matrix [Q] (plane stress)
Q11 = E1 / (1 - nu12 * nu21);
Q22 = E2 / (1 - nu12 * nu21);
Q12 = nu12 * E2 / (1 - nu12 * nu21);
Q66 = G12;
Q = [Q11, Q12, 0;
     Q12, Q22, 0;
     0,   0,   Q66];

% Ply interface coordinates through the thickness
z = linspace(-h/2, h/2, n_layers + 1);

% Laminate stiffness matrices
A = zeros(3, 3);
B = zeros(3, 3);
D = zeros(3, 3);

% Assemble [A], [B], [D]
for k = 1:n_layers
    % Transform ply stiffness to laminate axes
    th = deg2rad(theta(k));
    c = cos(th);
    s = sin(th);
    T = [c^2, s^2, 2*c*s;
         s^2, c^2, -2*c*s;
         -c*s, c*s, c^2 - s^2];
    T_inv = inv(T);
    Qbar = T_inv * Q * (T_inv)';

    % Thickness integration
    A = A + Qbar * (z(k+1) - z(k));
    B = B + 0.5 * Qbar * (z(k+1)^2 - z(k)^2);
    D = D + (1/3) * Qbar * (z(k+1)^3 - z(k)^3);
end

% Print stiffness matrices ([B] = 0 for a symmetric layup - quick sanity check)
fprintf('Matrix [A] (N/m):\n');
disp(A);
fprintf('Matrix [B] (N):\n');
disp(B);
fprintf('Matrix [D] (N*m):\n');
disp(D);

% Deflection by the Navier method
m_max = 5;        % series terms in x
n_max = 5;        % series terms in y
w_mn = zeros(m_max, n_max);  % deflection coefficients
x = a/2;          % evaluation point (plate center)
y = b/2;
w = 0;            % deflection

for m = 1:2:m_max            % odd m only (uniform load)
    for n = 1:2:n_max        % odd n only
        % Load coefficient q_mn
        q_mn = 16 * q0 / (pi^2 * m * n);

        % Stiffness denominator Delta_mn
        Delta_mn = D(1,1) * (m * pi / a)^4 + ...
                   2 * (D(1,2) + 2 * D(3,3)) * (m * pi / a)^2 * (n * pi / b)^2 + ...
                   D(2,2) * (n * pi / b)^4;

        % Deflection coefficient
        w_mn(m,n) = q_mn / Delta_mn;

        % Series summation
        w = w + w_mn(m,n) * sin(m * pi * x / a) * sin(n * pi * y / b);
    end
end

fprintf('Maximum center deflection (mm): %.4f\n', w * 1000);

% Curvatures at the center
kappa_x = 0;
kappa_y = 0;
kappa_xy = 0;
for m = 1:2:m_max
    for n = 1:2:n_max
        kappa_x = kappa_x + w_mn(m,n) * (m * pi / a)^2 * sin(m * pi * x / a) * sin(n * pi * y / b);
        kappa_y = kappa_y + w_mn(m,n) * (n * pi / b)^2 * sin(m * pi * x / a) * sin(n * pi * y / b);
        kappa_xy = kappa_xy + 2 * w_mn(m,n) * (m * pi / a) * (n * pi / b) * cos(m * pi * x / a) * cos(n * pi * y / b);
    end
end
kappa_x = -kappa_x;   % sign convention: kappa = -d2w/dx2
kappa_y = -kappa_y;
kappa = [kappa_x; kappa_y; kappa_xy];

fprintf('Curvatures at the center (1/m):\n');
fprintf('kappa_x = %.4e\n', kappa_x);
fprintf('kappa_y = %.4e\n', kappa_y);
fprintf('kappa_xy = %.4e\n', kappa_xy);

% Bending moments
M = D * kappa;
fprintf('Moments (N*m/m):\n');
fprintf('M_x = %.2f\n', M(1));
fprintf('M_y = %.2f\n', M(2));
fprintf('M_xy = %.2f\n', M(3));

% Ply-by-ply stresses
fprintf('\nPly stresses (MPa):\n');
for k = 1:n_layers
    fprintf('Ply %d (theta = %d deg):\n', k, theta(k));
    th = deg2rad(theta(k));
    c = cos(th);
    s = sin(th);
    T = [c^2, s^2, 2*c*s;
         s^2, c^2, -2*c*s;
         -c*s, c*s, c^2 - s^2];
    T_inv = inv(T);
    Qbar = T_inv * Q * (T_inv)';

    % Stress at the ply top surface (z = z(k+1))
    fprintf('  Top surface (z = %.4f m):\n', z(k+1));
    sigma_global = Qbar * kappa * z(k+1);
    sigma_local = T * sigma_global;   % global(xy) -> material(12): {sigma}_12 = [T]{sigma}_xy
    fprintf('    sigma_1 = %.2f MPa\n', sigma_local(1) / 1e6);
    fprintf('    sigma_2 = %.2f MPa\n', sigma_local(2) / 1e6);
    fprintf('    tau_12 = %.2f MPa\n', sigma_local(3) / 1e6);

    % Stress at the ply mid-surface
    fprintf('  Mid surface (z = %.4f m):\n', (z(k)+z(k+1))/2);
    sigma_global = Qbar * kappa * (z(k)+z(k+1))/2;
    sigma_local = T * sigma_global;
    fprintf('    sigma_1 = %.2f MPa\n', sigma_local(1) / 1e6);
    fprintf('    sigma_2 = %.2f MPa\n', sigma_local(2) / 1e6);
    fprintf('    tau_12 = %.2f MPa\n', sigma_local(3) / 1e6);

    % Stress at the ply bottom surface (z = z(k))
    fprintf('  Bottom surface (z = %.4f m):\n', z(k));
    sigma_global = Qbar * kappa * z(k);
    sigma_local = T * sigma_global;
    fprintf('    sigma_1 = %.2f MPa\n', sigma_local(1) / 1e6);
    fprintf('    sigma_2 = %.2f MPa\n', sigma_local(2) / 1e6);
    fprintf('    tau_12 = %.2f MPa\n', sigma_local(3) / 1e6);
end
