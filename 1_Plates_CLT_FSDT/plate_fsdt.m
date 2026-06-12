% Bending of a laminated composite plate - First-order Shear Deformation
% Theory (FSDT, Mindlin-Reissner).
%
% Same [0/90/0] plate as plate_clt.m, but with transverse shear flexibility.
% For each Navier harmonic (m,n) the COUPLED 3x3 system for (W, phi_x, phi_y)
% is solved - deflection and rotations cannot be decoupled in FSDT.
%
% Sanity check: in the thin-plate limit (large shear stiffness) the results
% converge to CLT; for this 6 mm plate the two theories agree closely.
clear all;
clc;

% Plate parameters
a = 1;          % plate length, m
b = 0.5;        % plate width, m
h = 0.006;      % total laminate thickness, m
h_k = 0.002;    % ply thickness, m
theta = [0, 90, 0];        % ply orientation angles, deg
n_layers = length(theta);  % number of plies
q0 = 10e3;      % uniform pressure, N/m^2
k_shear = 5/6;  % shear correction factor (rectangular section)

% Material properties
E1 = 140e9; E2 = 10e9; G12 = 5e9; nu12 = 0.3;
nu21 = nu12 * E2 / E1;
G13 = 5e9; G23 = 5e9;

% Ply stiffness matrix [Q] (plane stress)
Q11 = E1 / (1 - nu12 * nu21);
Q22 = E2 / (1 - nu12 * nu21);
Q12 = nu12 * E2 / (1 - nu12 * nu21);
Q66 = G12;
Q = [Q11, Q12, 0;
     Q12, Q22, 0;
     0,   0,   Q66];

% Transverse shear stiffness [Q_s]
Qs = [G13, 0;
      0,   G23];

% Ply interface coordinates
z = linspace(-h/2, h/2, n_layers + 1);

% Laminate stiffness matrices
A = zeros(3, 3);
B = zeros(3, 3);
D = zeros(3, 3);
A_s = zeros(2, 2);

for k = 1:n_layers
    th = deg2rad(theta(k));
    c = cos(th); s = sin(th);
    T = [c^2, s^2, 2*c*s;
         s^2, c^2, -2*c*s;
         -c*s, c*s, c^2 - s^2];
    T_inv = inv(T);
    Qbar = T_inv * Q * (T_inv)';

    % Transverse shear transforms with the plain 2x2 rotation matrix
    T_s = [ c, s;
           -s, c];
    Qs_bar = T_s * Qs * T_s';

    A   = A   + Qbar * (z(k+1) - z(k));
    B   = B   + 0.5 * Qbar * (z(k+1)^2 - z(k)^2);
    D   = D   + (1/3) * Qbar * (z(k+1)^3 - z(k)^3);
    A_s = A_s + k_shear * Qs_bar * (z(k+1) - z(k));
end

% Print stiffness matrices
fprintf('Matrix [A] (N/m):\n');  disp(A);
fprintf('Matrix [B] (N):\n');    disp(B);
fprintf('Matrix [D] (N*m):\n');  disp(D);
fprintf('Transverse shear stiffness [A_s] (N/m):\n'); disp(A_s);

% Shorthand
D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);
A55 = A_s(1,1);   % transverse shear xz
A44 = A_s(2,2);   % transverse shear yz

% --- Navier solution (FSDT): coupled 3x3 system per harmonic (m,n) ---
% w     = sum  W_mn  sin(alpha x) sin(beta y)
% phi_x = sum  X_mn  cos(alpha x) sin(beta y)
% phi_y = sum  Y_mn  sin(alpha x) cos(beta y)
m_max = 5; n_max = 5;
x = a/2; y = b/2;                 % evaluation point (center)
w = 0; phi_x = 0; phi_y = 0;
kappa_x = 0; kappa_y = 0; kappa_xy = 0;

for m = 1:2:m_max
    for n = 1:2:n_max
        al = m*pi/a; be = n*pi/b;
        q_mn = 16 * q0 / (pi^2 * m * n);

        s11 = A55*al^2 + A44*be^2;
        s12 = A55*al;
        s13 = A44*be;
        s22 = D11*al^2 + D66*be^2 + A55;
        s23 = (D12 + D66)*al*be;
        s33 = D66*al^2 + D22*be^2 + A44;

        Smat = [s11 s12 s13;
                s12 s22 s23;
                s13 s23 s33];
        sol = Smat \ [q_mn; 0; 0];
        W = sol(1); X = sol(2); Y = sol(3);

        sa = sin(m*pi/2); sb = sin(n*pi/2);
        ca = cos(m*pi/2); cb = cos(n*pi/2);

        w     = w     + W * sa * sb;
        phi_x = phi_x + X * ca * sb;
        phi_y = phi_y + Y * sa * cb;

        % Curvatures (sign convention consistent with CLT:
        % bottom fibers in tension under downward load)
        kappa_x  = kappa_x  + X * al * sa * sb;
        kappa_y  = kappa_y  + Y * be * sa * sb;
        kappa_xy = kappa_xy - (X * be + Y * al) * ca * cb;
    end
end

fprintf('Maximum center deflection (mm): %.4f\n', w * 1000);
fprintf('Rotation phi_x at center (rad): %.4e\n', phi_x);
fprintf('Rotation phi_y at center (rad): %.4e\n', phi_y);

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
    c = cos(th); s = sin(th);
    T = [c^2, s^2, 2*c*s;
         s^2, c^2, -2*c*s;
         -c*s, c*s, c^2 - s^2];
    T_inv = inv(T);
    Qbar = T_inv * Q * (T_inv)';

    fprintf('  Top surface (z = %.4f m):\n', z(k+1));
    sigma_global = Qbar * kappa * z(k+1);
    sigma_local = T * sigma_global;   % global(xy) -> material(12)
    fprintf('    sigma_1 = %.2f MPa\n', sigma_local(1) / 1e6);
    fprintf('    sigma_2 = %.2f MPa\n', sigma_local(2) / 1e6);
    fprintf('    tau_12 = %.2f MPa\n', sigma_local(3) / 1e6);

    fprintf('  Bottom surface (z = %.4f m):\n', z(k));
    sigma_global = Qbar * kappa * z(k);
    sigma_local = T * sigma_global;
    fprintf('    sigma_1 = %.2f MPa\n', sigma_local(1) / 1e6);
    fprintf('    sigma_2 = %.2f MPa\n', sigma_local(2) / 1e6);
    fprintf('    tau_12 = %.2f MPa\n', sigma_local(3) / 1e6);
end
