% Matrix crack density evolution - Discrete Damage Mechanics (DDM).
%
% Model in the spirit of Barbero & Cortes, Composites Part B 41 (2010):
% transverse tension opens matrix cracks (mode I, fed by eps2 through the
% degradation of Q22); in-plane shear contributes mode II (gamma12 through
% Q66). For each applied strain the crack density lambda follows from the
% damage activation condition g(lambda) = 0, solved by bisection (g is
% monotonically decreasing in lambda, so the root is unique and the
% iteration cannot fail).
%
% Two physical ingredients beyond the textbook minimum:
%
%   1. The damage law SATURATES: D2(inf) = D2max < 1. At the characteristic
%      damage state (CDS) matrix cracks unload the ply only partially -
%      the rest of the transverse stiffness is carried around the cracks
%      by shear transfer from the neighboring plies (shear lag). Without
%      saturation the ply stiffness would spuriously drop to zero.
%
%   2. Fracture toughness is STATISTICAL (Weibull, modulus m_w): cracks
%      nucleate at matrix defects with scattered strength. A single
%      deterministic G_IC produces a sharp kink of the stiffness curve
%      exactly at the threshold and a convex ~1/eps^2 decay - the smooth
%      upper shoulder seen in experiments (e.g. Varna & Berglund) is
%      mathematically impossible. Averaging over deterministic Weibull
%      quantiles (no random numbers, fully reproducible) gives the
%      experimentally observed descending sigmoid:
%      plateau 1 -> smooth knee -> inflection -> plateau 1-D2max.
%
% C_err and the law parameters (D2max, k2, D6max, k6, m_w) are calibration
% constants; set N_w = 1 to recover the deterministic model.
clc; clear; close all;

%% === CONSTANTS AND PARAMETERS ===
E1  = 140e9;   % Pa
E2  = 10e9;    % Pa
G12 = 5e9;     % Pa
nu12 = 0.3;
nu21 = nu12 * E2 / E1;
den  = 1 - nu12 * nu21;          % plane stress

G_IC  = 200;     % J/m^2, mode I fracture toughness (MEAN value)
G_IIC = 600;     % J/m^2, mode II fracture toughness (MEAN value)
r = G_IC / G_IIC;   % mode mixity ratio (independent of the scatter, see below)

% Strength statistics (Weibull): site i has toughness G_IC*xi_i, G_IIC*xi_i
% (a defect weakens both modes equally => r is preserved).
% Deterministic quantiles: xi_i = scale*(-ln(1-F_i))^(1/m_w), F_i=(i-0.5)/N_w;
% scale normalizes the mean, E[xi] = 1. N_w = 1 -> deterministic model.
m_w = 4;         % typical Weibull modulus for matrix defects
N_w = 40;
if N_w == 1
    xi_w = 1;
else
    Fq   = ((1:N_w) - 0.5) / N_w;
    xi_w = (-log(1 - Fq)).^(1/m_w) / gamma(1 + 1/m_w);
end

h = 0.000125;    % m, ply thickness
C_err = 10 * h;  % m, energy length scale of the cracking cell (calibration)

gamma12 = 1e-3;  % fixed shear strain (mode II contribution)

strain_steps = linspace(0, 0.015, 150);   % applied transverse strain eps2

%% === DAMAGE LAWS AND ENERGY RELEASE RATES ===
% Damage parameters: 0 at lambda = 0, SATURATION at Dmax < 1 as lambda -> inf
% (CDS: cracks remove only part of the ply stiffness, the rest is shear lag).
% The product Dmax*k sets the slope dD/dlambda at zero and thereby the
% cracking threshold.
D2max = 0.7;  k2 = 2.0 / D2max;        % saturation and rate, mode I
D6max = 0.8;  k6 = 1.2 / D6max;        % same for mode II
D2  = @(lam) D2max * (1 - exp(-k2 * lam));
D6  = @(lam) D6max * (1 - exp(-k6 * lam));
dD2 = @(lam) D2max * k2 * exp(-k2 * lam);    % dD2/dlambda
dD6 = @(lam) D6max * k6 * exp(-k6 * lam);    % dD6/dlambda

% Energy release rates: G = -(C_err/2) * eps^2 * dQ/dlambda  (>= 0)
% dQ22/dlambda = -E2*dD2/den  ->  G_I  = (C_err/2)*eps2^2*E2*dD2/den
% dQ66/dlambda = -G12*dD6     ->  G_II = (C_err/2)*gam^2 *G12*dD6
GI  = @(lam, e2)  0.5 * C_err * e2.^2  * (E2  * dD2(lam) / den);
GII = @(lam, gg)  0.5 * C_err * gg.^2  * (G12 * dD6(lam));

% Damage activation function for a site with toughness fraction xw
% (g <= 0: no crack growth; g = 0: equilibrium)
gfun = @(lam, e2, gg, xw) (1 - r) * sqrt(GI(lam, e2) / (G_IC * xw)) ...
       + r * (GI(lam, e2) / (G_IC * xw) + GII(lam, gg) / (G_IIC * xw)) - 1;

%% === SOLUTION: g(lambda) = 0 by bisection at each strain level,
%%     averaged over the Weibull quantiles ===
N = length(strain_steps);
lambda_history = zeros(1, N);   % mean crack density over the sites
GI_history     = zeros(1, N);
GII_history    = zeros(1, N);
D2_history     = zeros(1, N);   % mean damage (this is what drives stiffness)

lams = zeros(1, numel(xi_w));
for i = 1:N
    eps2 = strain_steps(i);
    for k = 1:numel(xi_w)
        xw = xi_w(k);
        if gfun(0, eps2, gamma12, xw) <= 0
            lams(k) = 0;              % threshold of this site not reached
        else
            % g(lambda) decreases monotonically: find an upper bound with g<0
            lo = 0; hi = 1;
            while gfun(hi, eps2, gamma12, xw) > 0 && hi < 1e4
                hi = hi * 2;
            end
            for it = 1:80             % bisection (hi*2^-80 << problem scale)
                mid = 0.5 * (lo + hi);
                if gfun(mid, eps2, gamma12, xw) > 0
                    lo = mid;
                else
                    hi = mid;
                end
            end
            lams(k) = 0.5 * (lo + hi);
        end
    end
    lambda_history(i) = mean(lams);
    D2_history(i)     = mean(D2(lams));          % average the DAMAGE
    GI_history(i)     = GI(lambda_history(i), eps2);     % diagnostics at
    GII_history(i)    = GII(lambda_history(i), gamma12); % the mean lambda
end

%% === PLOTTING ===
figure;
subplot(3,1,1);
plot(strain_steps*100, lambda_history, 'b', 'LineWidth', 2);
xlabel('\epsilon_2 (%)'); ylabel('\lambda (crack density)');
title('Crack Density Evolution'); grid on;

subplot(3,1,2);
plot(strain_steps*100, GI_history, 'r', 'LineWidth', 2); hold on;
plot(strain_steps*100, GII_history, 'g', 'LineWidth', 2);
yline(G_IC,  'r--', 'G_{IC}');
yline(G_IIC, 'g--', 'G_{IIC}');
xlabel('\epsilon_2 (%)'); ylabel('G_I, G_{II} (J/m^2)');
legend('G_I', 'G_{II}', 'Location', 'best');
title('Energy Release Rates'); grid on;

subplot(3,1,3);
plot(strain_steps*100, 1 - D2_history, 'k', 'LineWidth', 2);
yline(1 - D2max, 'k--', 'CDS plateau: 1-D_{2max}');
xlabel('\epsilon_2 (%)'); ylabel('E_2^{eff}/E_2 = 1-D_2');
title('Transverse Stiffness Reduction (saturating sigmoid)');
grid on; ylim([0 1.05]);

dR = gradient(1 - D2_history, strain_steps);
[~, i_inf] = min(dR);                 % sigmoid inflection (steepest descent)
fprintf('First cracks (weakest sites): eps2 = %.2f %%\n', ...
        strain_steps(find(lambda_history > 1e-3, 1)) * 100);
fprintf('Stiffness curve inflection: eps2 = %.2f %%\n', strain_steps(i_inf)*100);
fprintf('Maximum crack density lambda = %.3f at eps2 = %.2f %%\n', ...
        lambda_history(end), strain_steps(end)*100);
fprintf('Residual stiffness E2eff/E2 = %.3f (CDS plateau: %.2f)\n', ...
        1 - D2_history(end), 1 - D2max);
