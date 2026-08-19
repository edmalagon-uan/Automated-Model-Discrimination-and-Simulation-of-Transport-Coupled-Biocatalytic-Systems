%% immobilized_catalase_reusability.m
%% REUSABILITY KINETIC SIMULATOR: IMMOBILIZED CATALASE (Kaushal et al., 2018)
% Simulates the cycle-by-cycle activity loss of catalase immobilized on
% chitosan (CS) and on a chitosan-bentonite complex (CS-B), using a
% first-order activity-decay model calibrated with the activity retentions
% reported in the reference article (50% for CS and 70% for CS-B after
% 20 reuse cycles).
%
% Produces the figure used as Figure 4 of the SSICE manuscript.
%
% Requirements: MATLAB R2018b or later (yline); no toolboxes required
%
% Reference: Kaushal, J., Singh, G., & Arya, S. K. (2018). Immobilization
% of catalase onto chitosan and chitosan-bentonite complex: A comparative
% study. Biotechnology Reports, 18, e00258.
% https://doi.org/10.1016/j.btre.2018.e00258
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

clear; clc; close all;
fprintf('=========================================================================\n');
fprintf(' IMMOBILIZED CATALASE KINETICS AND REUSABILITY SIMULATOR\n');
fprintf('=========================================================================\n\n');

% =========================================================================
% 1. KINETIC PARAMETERS FROM THE LITERATURE (Kaushal et al., 2018)
% =========================================================================
% Initial substrate concentration for the assay (H2O2), mM
S_assay = 50.0;

% Michaelis-Menten constants (Km) reported in the paper (mM)
Km_free = 12.5;
Km_CS   = 25.0; % Catalase on pure chitosan
Km_CSB  = 20.0; % Catalase on chitosan-bentonite

% Theoretical initial maximum velocity (relative units based on the paper)
Vmax_initial = 33000;

% Number of reuse cycles evaluated
n_cycles = 20;
cycles = 1:n_cycles;

% =========================================================================
% 2. CALCULATION OF THE DEACTIVATION RATE (MATHEMATICAL MODELING)
% =========================================================================
% The article reports that after 20 cycles:
% - Pure chitosan (CS) retains 50% (0.50)
% - Chitosan-bentonite (CS-B) retains 70% (0.70)
% We use the exponential activity-decay model: Activity = Exp(-k * cycle)
k_deact_CS  = -log(0.50) / n_cycles;
k_deact_CSB = -log(0.70) / n_cycles;

fprintf('>> Deactivation rates calculated by the model:\n');
fprintf('   Pure chitosan support (CS)       | k = %.4f per cycle\n', k_deact_CS);
fprintf('   Chitosan-bentonite support (CSB) | k = %.4f per cycle\n\n', k_deact_CSB);

% =========================================================================
% 3. CYCLE-BY-CYCLE KINETIC VELOCITY SIMULATION
% =========================================================================
v_CS  = zeros(1, n_cycles);
v_CSB = zeros(1, n_cycles);

for c = 1:n_cycles
    % 1. Compute the remaining Vmax in this specific cycle due to wear
    Vmax_CS_current  = Vmax_initial * exp(-k_deact_CS * c);
    Vmax_CSB_current = Vmax_initial * exp(-k_deact_CSB * c);

    % 2. Apply the kinetic equation corresponding to each support
    % Note how each one uses the Km altered by the diffusional support
    v_CS(c)  = (Vmax_CS_current  * S_assay) / (Km_CS  + S_assay);
    v_CSB(c) = (Vmax_CSB_current * S_assay) / (Km_CSB + S_assay);
end

% Free-enzyme velocity as reference control (Cycle 1)
v_free_control = (Vmax_initial * S_assay) / (Km_free + S_assay);

% =========================================================================
% 4. PROFESSIONAL GRAPHICAL DISPLAY
% =========================================================================
figure('Name', 'Catalase Stability and Reusability', 'Color', 'w', 'Position', [200, 200, 700, 450]);

% Plot the behavior of the pure chitosan support
plot(cycles, v_CS, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'Catalase on Chitosan (CS) [50% retention]');
hold on;

% Plot the behavior of the chitosan-bentonite support
plot(cycles, v_CSB, 'b^--', 'LineWidth', 2, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'Catalase on Chitosan-Bentonite (CS-B) [70% retention]');

% Horizontal reference line of the free enzyme (lost after 1 cycle)
yline(v_free_control, 'k:', 'LineWidth', 1.5, ...
    'DisplayName', 'Free Enzyme (initial control, no support)');

% Aesthetic adjustments of the plot
xlabel('Reuse Cycle Number', 'FontSize', 11);
ylabel('Estimated Reaction Velocity (v)', 'FontSize', 11);
title('Simulation of Catalase Activity Loss upon Reuse', 'FontSize', 12);
xticks(1:2:n_cycles); % Show cycles in steps of 2 for visual clarity
legend('Location', 'NorthEast', 'FontSize', 10);
grid on;

% =========================================================================
% 5. PREDICTIVE REPORT OF INDUSTRIAL PERFORMANCE
% =========================================================================
fprintf('--- PREDICTIVE SUMMARY OF THE IMMOBILIZED SYSTEM ---\n');
fprintf('Estimated initial velocity (Cycle 1):\n');
fprintf('  - Free enzyme: %.2f\n', v_free_control);
fprintf('  - On chitosan (CS): %.2f (affected by high Km)\n', v_CS(1));
fprintf('  - On matrix (CS-B): %.2f (best initial performance)\n\n', v_CSB(1));
fprintf('Long-term loss prediction (Cycle 20):\n');
fprintf('  - Velocity on pure chitosan drops to: %.2f\n', v_CS(end));
fprintf('  - Velocity on chitosan-bentonite drops to: %.2f\n', v_CSB(end));
fprintf('>> MODEL CONCLUSION: The chitosan-bentonite hybrid composite protects the\n');
fprintf('   catalase structure better against mechanical deactivation in continuous cycles.\n');
