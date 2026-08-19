%% batch_reactor_product_inhibition.m
%% REACTOR ENGINEERING INTEGRATION MODULE
% Simulation of a progress curve in a batch reactor with uncompetitive
% product inhibition.
% Based on the kinetic modeling for maltose hydrolysis (Bas et al., 2007).
%
% Produces the figure used as Figure 3 of the SSICE manuscript.
%
% Requirements: MATLAB (base installation; no toolboxes required)
%
% Reference: Bas, D., Duduk, F. C., & Boyaci, I. H. (2007). Modeling and
% optimization IV: Investigation of reaction kinetics and kinetic
% parameters using artificial neural networks. International Journal of
% Food Engineering, 3(1).
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

clear; clc; close all;

% 1. Kinetic parameters obtained from the article's model
Vmax = 1.48;   % Maximum velocity (umol/min/mg enzyme)
Km   = 1.91;   % Michaelis-Menten constant (mM)
Ki   = 71.42;  % Uncompetitive product inhibition constant (mM)

% 2. Initial conditions of the batch reactor
S0      = 16.0; % Initial substrate concentration [Maltose] (mM)
P0      = 0.0;  % Initial product concentration [Glucose] (mM)
t_final = 80;   % Total simulation operation time (minutes)
dt      = 0.1;  % Time-step size for the numerical method (min)

% 3. Initialization of vectors to store the time simulation
time_sim = 0:dt:t_final;
S_sim = zeros(size(time_sim));
P_sim = zeros(size(time_sim));

% Assignment of values at time t = 0
S_sim(1) = S0;
P_sim(1) = P0;

% 4. Solution of the differential equation (finite differences / Euler method)
for i = 1:(length(time_sim)-1)
    S_current = S_sim(i);
    P_current = P_sim(i);

    % Mathematical equation of the uncompetitive product inhibition model
    % The article assumes the stoichiometric relation P = P0 + 2*(S0 - S)
    v_reaction = (Vmax * S_current) / (Km + S_current * (1 + P_current / Ki));

    % Step-by-step numerical approximation: S(t+1) = S(t) - v * dt
    S_sim(i+1) = S_current - (v_reaction * dt);

    % Physical safety constraint: chemical substrate cannot be below zero
    if S_sim(i+1) < 0
        S_sim(i+1) = 0;
    end

    % Stoichiometric update of the product generated (1 Maltose -> 2 Glucose)
    P_sim(i+1) = P0 + 2 * (S0 - S_sim(i+1));
end

% 5. Graphical display of the reactor kinetics
figure('Name', 'Reactor Simulation - Bas et al. 2007 Case', 'Color', 'w');
plot(time_sim, S_sim, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Substrate [S] (Maltose)'); hold on;
plot(time_sim, P_sim, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Product [P] (Glucose)');

% Aesthetic configuration of the plot
xlabel('Reaction Time (minutes)', 'FontSize', 11);
ylabel('Concentration in the Reactor (mM)', 'FontSize', 11);
title('Dynamic Simulation of the Batch Reactor (Neural Network Constants)', 'FontSize', 12);
legend('Location', 'East', 'FontSize', 10);
grid on;

% Print summary to the command window
fprintf('=== REACTOR SIMULATION COMPLETED ===\n');
fprintf('Evaluated mechanism: Uncompetitive Product Inhibition\n');
fprintf('Final substrate concentration: %.4f mM at %d minutes.\n', S_sim(end), t_final);
fprintf('Final product concentration: %.4f mM.\n', P_sim(end));
