% Sectional load analysis of a rocket body: axial force N, shear force Q
% and bending moment M diagrams along the length, plus pressure-vessel
% stresses in the propellant tanks (cylinder and bottom, with automatic
% selection of the lighter bottom type).

% Rocket dimensions
L = 10; % length, m
n_sections = 100; % number of calculation sections
x = linspace(0, L, n_sections); % coordinates of sections, m
D = ones(size(x)) * 1; % diameters at sections, m (assuming constant)

% Tank types: 1 - integrated, 0 - suspended
tank_type_fuel = 1;
tank_type_oxidizer = 1;

% Tank coordinates
if tank_type_fuel == 1
    x_start_fuel = 1; % start of fuel tank, m
    x_end_fuel = 3; % end of fuel tank, m
else
    x_attach_fuel = 1.5; % attachment point of suspended fuel tank, m
end

if tank_type_oxidizer == 1
    x_start_oxidizer = 4; % start of oxidizer tank, m
    x_end_oxidizer = 6; % end of oxidizer tank, m
else
    x_attach_oxidizer = 4.5; % attachment point of suspended oxidizer tank, m
end

% Tank wall thickness, m
delta_fuel = 0.005;
delta_oxidizer = 0.005;

% Pressurization pressure, Pa
p_fuel = 1e5;
p_oxidizer = 1e5;

% Cargo
x_cm_cargo = [2, 5]; % coordinates of cargo centers of mass, m
m_cargo = [100, 200]; % masses of cargo, kg
x_attach_cargo = [2, 5]; % attachment points of cargo, m

% Engine
x_cm_engine = 9; % coordinate of engine center of mass, m
m_engine = 500; % mass of engine, kg
x_thrust = 10; % point of thrust application, m

% Thrust
T = 10000; % thrust magnitude, N
theta = 0; % thrust angle, rad (0 for axial thrust)

% Flight parameters
h = 10000; % altitude, m
gamma = 0; % pitch angle, rad
m_total = 10000; % total flight mass, kg
m_fuel = 2000; % mass of fuel, kg
m_oxidizer = 3000; % mass of oxidizer, kg

g = 9.81; % gravity, m/s^2

% Calculate structural mass
m_structure = m_total - m_fuel - m_oxidizer - sum(m_cargo) - m_engine;

% Stress calculation functions are defined at the END of this file
% (MATLAB requires local functions in a script to appear after all script code).

% Select bottom types for integrated tanks
if tank_type_fuel == 1
    R_fuel = D(find(x >= x_start_fuel, 1)) / 2;
    sigma_spherical_fuel = calc_spherical_stress(p_fuel, R_fuel, delta_fuel);
    sigma_flat_fuel = calc_flat_stress(p_fuel, R_fuel, delta_fuel);
    if sigma_spherical_fuel < sigma_flat_fuel
        bottom_type_fuel = 'spherical';
    else
        bottom_type_fuel = 'flat';
    end
end

if tank_type_oxidizer == 1
    R_oxidizer = D(find(x >= x_start_oxidizer, 1)) / 2;
    sigma_spherical_oxidizer = calc_spherical_stress(p_oxidizer, R_oxidizer, delta_oxidizer);
    sigma_flat_oxidizer = calc_flat_stress(p_oxidizer, R_oxidizer, delta_oxidizer);
    if sigma_spherical_oxidizer < sigma_flat_oxidizer
        bottom_type_oxidizer = 'spherical';
    else
        bottom_type_oxidizer = 'flat';
    end
end

% Distributed load from structure
q_structure = (m_structure / L) * g * cos(gamma); % N/m

% Initialize total distributed load
q = q_structure * ones(size(x));

% Add distributed loads from integrated tanks
if tank_type_fuel == 1
    idx_fuel = (x >= x_start_fuel & x <= x_end_fuel);
    L_fuel = x_end_fuel - x_start_fuel;
    q_fuel = (m_fuel / L_fuel) * g * cos(gamma);
    q(idx_fuel) = q(idx_fuel) + q_fuel;
end

if tank_type_oxidizer == 1
    idx_oxidizer = (x >= x_start_oxidizer & x <= x_end_oxidizer);
    L_oxidizer = x_end_oxidizer - x_start_oxidizer;
    q_oxidizer = (m_oxidizer / L_oxidizer) * g * cos(gamma);
    q(idx_oxidizer) = q(idx_oxidizer) + q_oxidizer;
end

% Concentrated loads from suspended tanks
P_suspended = zeros(size(x));
if tank_type_fuel == 0
    idx_attach_fuel = find(x >= x_attach_fuel, 1);
    P_suspended(idx_attach_fuel) = P_suspended(idx_attach_fuel) - m_fuel * g * cos(gamma);
end
if tank_type_oxidizer == 0
    idx_attach_oxidizer = find(x >= x_attach_oxidizer, 1);
    P_suspended(idx_attach_oxidizer) = P_suspended(idx_attach_oxidizer) - m_oxidizer * g * cos(gamma);
end

% Concentrated loads from cargo
P_cargo = zeros(size(x));
for i = 1:length(x_attach_cargo)
    idx = find(x >= x_attach_cargo(i), 1);
    P_cargo(idx) = P_cargo(idx) - m_cargo(i) * g * cos(gamma);
end

% Concentrated load from engine
P_engine = zeros(size(x));
idx_engine = find(x >= x_cm_engine, 1);
P_engine(idx_engine) = P_engine(idx_engine) - m_engine * g * cos(gamma);

% Thrust force (transverse component)
P_thrust = zeros(size(x));
idx_thrust = find(x >= x_thrust, 1);
P_thrust(idx_thrust) = P_thrust(idx_thrust) + T * sin(theta);

% Axial force due to pressure for integrated tanks with flat bottoms
N_pressure = zeros(size(x));
if tank_type_fuel == 1 && strcmp(bottom_type_fuel, 'flat')
    idx_fuel = (x > x_start_fuel & x < x_end_fuel);
    R_fuel = D(find(x >= x_start_fuel, 1)) / 2;
    N_pressure(idx_fuel) = - p_fuel * pi * R_fuel^2;
end
if tank_type_oxidizer == 1 && strcmp(bottom_type_oxidizer, 'flat')
    idx_oxidizer = (x > x_start_oxidizer & x < x_end_oxidizer);
    R_oxidizer = D(find(x >= x_start_oxidizer, 1)) / 2;
    N_pressure(idx_oxidizer) = - p_oxidizer * pi * R_oxidizer^2;
end

% Calculate internal forces and moments
N = zeros(size(x));
Q = zeros(size(x));
M = zeros(size(x));

for i = 2:n_sections
    dx = x(i) - x(i-1);

    % Axial force N
    N(i) = N(i-1);
    if x(i-1) < x_thrust && x(i) >= x_thrust
        N(i) = N(i) + T * cos(theta);
    end
    % Add pressure effect
    N(i) = N(i) + N_pressure(i);

    % Transverse force Q
    Q(i) = Q(i-1) - q(i-1) * dx - P_suspended(i) - P_cargo(i) - P_engine(i) - P_thrust(i);

    % Bending moment M
    M(i) = M(i-1) + Q(i-1) * dx - (P_suspended(i) + P_cargo(i) + P_engine(i) + P_thrust(i)) * (dx / 2);
end

% Plot results
figure;
subplot(3,1,1);
plot(x, N, 'b-', 'LineWidth', 1.5);
title('Axial Force N, N');
xlabel('x, m');
ylabel('N, N');
grid on;

subplot(3,1,2);
plot(x, Q, 'r-', 'LineWidth', 1.5);
title('Shear Force Q, N');
xlabel('x, m');
ylabel('Q, N');
grid on;

subplot(3,1,3);
plot(x, M, 'g-', 'LineWidth', 1.5);
title('Bending Moment M, N*m');
xlabel('x, m');
ylabel('M, N*m');
grid on;

% Output selected bottom types and stresses
if tank_type_fuel == 1
    fprintf('Fuel tank bottom type: %s\n', bottom_type_fuel);
    [sigma_theta_fuel, sigma_x_fuel] = calc_cylinder_stress(p_fuel, R_fuel, delta_fuel);
    fprintf('Fuel tank cylinder: sigma_theta = %.2f Pa, sigma_x = %.2f Pa\n', sigma_theta_fuel, sigma_x_fuel);
    if strcmp(bottom_type_fuel, 'spherical')
        sigma_bottom_fuel = calc_spherical_stress(p_fuel, R_fuel, delta_fuel);
    else
        sigma_bottom_fuel = calc_flat_stress(p_fuel, R_fuel, delta_fuel);
    end
    fprintf('Fuel tank bottom: sigma = %.2f Pa\n', sigma_bottom_fuel);
end

if tank_type_oxidizer == 1
    fprintf('Oxidizer tank bottom type: %s\n', bottom_type_oxidizer);
    [sigma_theta_oxidizer, sigma_x_oxidizer] = calc_cylinder_stress(p_oxidizer, R_oxidizer, delta_oxidizer);
    fprintf('Oxidizer tank cylinder: sigma_theta = %.2f Pa, sigma_x = %.2f Pa\n', sigma_theta_oxidizer, sigma_x_oxidizer);
    if strcmp(bottom_type_oxidizer, 'spherical')
        sigma_bottom_oxidizer = calc_spherical_stress(p_oxidizer, R_oxidizer, delta_oxidizer);
    else
        sigma_bottom_oxidizer = calc_flat_stress(p_oxidizer, R_oxidizer, delta_oxidizer);
    end
    fprintf('Oxidizer tank bottom: sigma = %.2f Pa\n', sigma_bottom_oxidizer);
end

% ===== Local functions (must appear at the end of the script) =====
function sigma = calc_spherical_stress(p, R, delta)
    sigma = (p * R) / (2 * delta);
end

function sigma = calc_flat_stress(p, R, delta)
    sigma = (p * R^2) / (2 * delta); % simplified
end

function [sigma_theta, sigma_x] = calc_cylinder_stress(p, R, delta)
    sigma_theta = (p * R) / delta;
    sigma_x = (p * R) / (2 * delta);
end
