% Three-layer (sandwich) cylindrical shell with closed ends under internal
% pressure. The core shear stiffness enters through Gc.
%
% Solution: membrane state + edge-effect correction (clamped ends), the
% classic exponentially decaying boundary-layer solution. Membrane values
% follow the pressure-vessel formulas (sigma_theta = pR/(h1+h2),
% sigma_x = pR/2(h1+h2)) and serve as a built-in sanity check.
clear all; close all; clc;

% Parameters
R = 1;              % radius, m
L = 1;              % length, m
h1 = 0.01;          % outer face sheet thickness, m
h2 = 0.01;          % inner face sheet thickness, m
hc = 0.05;          % core thickness, m
E = 2e11;           % face sheet Young's modulus, Pa
nu = 0.3;           % Poisson's ratio
Gc = 1e8;           % core shear modulus, Pa
p = 2e6;            % internal pressure, Pa

% Stiffnesses
C = E * (h1 + h2) / (1 - nu^2);  % membrane stiffness
D = E * (h1^3 + h2^3) / (12 * (1 - nu^2)) + Gc * hc * (h1 + h2 + hc)^2 / 4;  % bending stiffness
beta = (D * R^2 / C)^(1/4);
alpha = 1 / (sqrt(2) * beta);

% Axial coordinate
x = linspace(0, L, 1000);

% Membrane solution
w_mem = p * R^2 / (E * (h1 + h2)) * (1 - nu/2);
u_mem = zeros(size(x));
sigma_x_mem = p * R / (2 * (h1 + h2));
sigma_theta_mem = p * R / (h1 + h2);
tau_xz_mem = zeros(size(x));
gamma_xz_mem = tau_xz_mem / Gc;

% Diagnostics (pressure-vessel check values)
disp(['Membrane sigma_theta = ', num2str(sigma_theta_mem / 1e6), ' MPa']);
disp(['Membrane sigma_x = ', num2str(sigma_x_mem / 1e6), ' MPa']);
disp(['Membrane w = ', num2str(w_mem * 1e3), ' mm']);

% Full solution (membrane + edge effect at both ends)
w0 = p * R^2 * (1 - nu^2) / (E * (h1 + h2));
w = w0 * (1 - (exp(-alpha * x) .* (cos(alpha * x) + sin(alpha * x)) + ...
               exp(-alpha * (L - x)) .* (cos(alpha * (L - x)) + sin(alpha * (L - x)))) / ...
               (1 + exp(-alpha * L)));
u = zeros(size(x)); % u(x) = 0 for u(0) = u(L) = 0

% Second derivative (bending)
d2w_dx2 = w0 * alpha^2 * (exp(-alpha * x) .* (cos(alpha * x) - sin(alpha * x)) - ...
          exp(-alpha * (L - x)) .* (cos(alpha * (L - x)) - sin(alpha * (L - x)))) / ...
          (1 + exp(-alpha * L));

% Face sheet stresses
z1 = hc/2 + h1/2;
sigma_x1 = p * R / (2 * (h1 + h2)) + E * z1 / (1 - nu^2) * d2w_dx2;
sigma_theta1 = p * R / (h1 + h2) * ones(size(x));

z2 = -hc/2 - h2/2;
sigma_x2 = p * R / (2 * (h1 + h2)) + E * z2 / (1 - nu^2) * d2w_dx2;
sigma_theta2 = p * R / (h1 + h2) * ones(size(x));

% Core shear stress
d3w_dx3 = -w0 * alpha^3 * (exp(-alpha * x) .* sin(alpha * x) + ...
           exp(-alpha * (L - x)) .* sin(alpha * (L - x))) / ...
           (1 + exp(-alpha * L));
tau_xz = D / hc * d3w_dx3;
tau_theta_z = zeros(size(x));
gamma_xz = tau_xz / Gc;

% Plots
figure('Position', [100, 100, 1200, 1000]);

% 1. Radial displacement
subplot(4, 2, 1);
h(1) = plot(x, w * 1e3, 'b-', 'LineWidth', 2, 'DisplayName', 'Full Solution');
hold on;
h(2) = plot(x, w_mem * ones(size(x)) * 1e3, 'r--', 'LineWidth', 2, 'DisplayName', 'Membrane Solution');
title('Radial Displacement w(x)');
xlabel('x (m)'); ylabel('w (mm)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [min(w * 1e3) * 0.9, max(w * 1e3) * 1.1]);

% 2. Axial displacement
subplot(4, 2, 2);
h(1) = plot(x, u * 1e3, 'b-', 'LineWidth', 2, 'DisplayName', 'Full Solution');
hold on;
h(2) = plot(x, u_mem * 1e3, 'r--', 'LineWidth', 2, 'DisplayName', 'Membrane Solution');
title('Axial Displacement u(x)');
xlabel('x (m)'); ylabel('u (mm)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [-0.1, 0.1]);

% 3. Axial stress (outer face)
subplot(4, 2, 3);
h(1) = plot(x, sigma_x1 / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\sigma_x^{(1)} (Full)');
hold on;
h(2) = plot(x, sigma_x_mem * ones(size(x)) / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', '\sigma_x^{(1)} (Membrane)');
title('Axial Stress \sigma_x^{(1)}');
xlabel('x (m)'); ylabel('\sigma_x^{(1)} (MPa)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [min(sigma_x1 / 1e6) * 1.1, max(sigma_x1 / 1e6) * 1.1]);

% 4. Circumferential stress (outer face)
subplot(4, 2, 4);
h(1) = plot(x, sigma_theta1 / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\sigma_\theta^{(1)} (Full)');
hold on;
h(2) = plot(x, sigma_theta_mem * ones(size(x)) / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', '\sigma_\theta^{(1)} (Membrane)');
title('Circumferential Stress \sigma_\theta^{(1)}');
xlabel('x (m)'); ylabel('\sigma_\theta^{(1)} (MPa)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [0, max(sigma_theta1 / 1e6) * 1.5]);

% 5. Axial stress (inner face)
subplot(4, 2, 5);
h(1) = plot(x, sigma_x2 / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\sigma_x^{(2)} (Full)');
hold on;
h(2) = plot(x, sigma_x_mem * ones(size(x)) / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', '\sigma_x^{(2)} (Membrane)');
title('Axial Stress \sigma_x^{(2)}');
xlabel('x (m)'); ylabel('\sigma_x^{(2)} (MPa)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [min(sigma_x2 / 1e6) * 1.1, max(sigma_x2 / 1e6) * 1.1]);

% 6. Circumferential stress (inner face)
subplot(4, 2, 6);
h(1) = plot(x, sigma_theta2 / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\sigma_\theta^{(2)} (Full)');
hold on;
h(2) = plot(x, sigma_theta_mem * ones(size(x)) / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', '\sigma_\theta^{(2)} (Membrane)');
title('Circumferential Stress \sigma_\theta^{(2)}');
xlabel('x (m)'); ylabel('\sigma_\theta^{(2)} (MPa)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [0, max(sigma_theta2 / 1e6) * 1.5]);

% 7. Core shear stress tau_xz
subplot(4, 2, 7);
h(1) = plot(x, tau_xz / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\tau_{xz} (Full)');
hold on;
h(2) = plot(x, tau_xz_mem / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', '\tau_{xz} (Membrane)');
title('Shear Stress \tau_{xz}');
xlabel('x (m)'); ylabel('\tau_{xz} (MPa)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [min(tau_xz / 1e6) * 1.1, max(tau_xz / 1e6) * 1.1]);

% 8. Shear stress tau_theta_z
subplot(4, 2, 8);
hh = plot(x, tau_theta_z / 1e6, 'b-', 'LineWidth', 2, 'DisplayName', '\tau_{\theta z} (Full)');
title('Shear Stress \tau_{\theta z}');
xlabel('x (m)'); ylabel('\tau_{\theta z} (MPa)');
legend(hh, {'Full Solution'}, 'Location', 'best');
grid on;

% Core shear strain
figure('Position', [100, 100, 600, 400]);
h(1) = plot(x, gamma_xz * 1e3, 'b-', 'LineWidth', 2, 'DisplayName', '\gamma_{xz} (Full)');
hold on;
h(2) = plot(x, gamma_xz_mem * 1e3, 'r--', 'LineWidth', 2, 'DisplayName', '\gamma_{xz} (Membrane)');
title('Shear Strain \gamma_{xz}');
xlabel('x (m)'); ylabel('\gamma_{xz} (mrad)');
legend(h, {'Full Solution', 'Membrane Solution'}, 'Location', 'best');
grid on;
set(gca, 'YLim', [min(gamma_xz * 1e3) * 1.1, max(gamma_xz * 1e3) * 1.1]);

sgtitle('Three-Layer Cylindrical Shell with Closed Ends');
