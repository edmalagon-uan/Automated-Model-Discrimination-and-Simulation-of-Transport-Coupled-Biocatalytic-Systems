%% ssice_integrated_kinetics.m
%% INTEGRATED PROFESSIONAL ENZYME KINETICS SUPER SCRIPT (SSICE)
% Developed for v0 processing, multi-model analysis and allosteric/
% inhibition diagnosis with automated model selection by AIC.
%
% Module I  : automated extraction of initial velocities (v0) from raw
%             time-vs-signal curves using a moving-window linearity
%             algorithm.
% Module II : macroscopic multi-model fitting. Without inhibitor, five
%             substrate-dependent models compete (Michaelis-Menten, Hill,
%             Haldane, MWC, substrate activation); with inhibitor, four
%             inhibition models compete (competitive, pure non-competitive,
%             uncompetitive, mixed). The best model is selected by the
%             Akaike Information Criterion (AIC).
%
% Produces the figures used as Figure 2 of the SSICE manuscript.
%
% Requirements: MATLAB R2016b or later
%               Optimization Toolbox (lsqcurvefit)
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

clear; clc; close all;
fprintf('=========================================================================\n');
fprintf(' STARTING INTEGRATED ENZYME KINETICS SUPER SCRIPT (SSICE)\n');
fprintf('=========================================================================\n\n');

% Global options for the mathematical optimization algorithms
options = optimset('Display','off','MaxFunEvals',20000,'MaxIter',20000);

%% =========================================================================
% MODULE I: AUTOMATED PROCESSING OF INITIAL VELOCITIES (v0)
% =========================================================================
disp('>> MODULE I: Processing raw time-vs-signal data...');
% We simulate time-course data for 3 experimental conditions as recorded
% by a plate reader or spectrophotometer.
% In practice, you can load this from an Excel file:
% time = matrix(:,1); signal = matrix(:,2);
common_time = (0:5:60)'; % Time every 5 seconds up to 1 minute

% Raw experimental curves (simulated signal/absorbance for 3 example assays)
assay_control_S1   = [0.002; 0.045; 0.088; 0.130; 0.171; 0.210; 0.245; 0.275; 0.300; 0.320; 0.335; 0.345; 0.350];
assay_control_S2   = [0.002; 0.085; 0.165; 0.240; 0.310; 0.370; 0.420; 0.460; 0.490; 0.515; 0.535; 0.550; 0.560];
assay_inhibitor_S2 = [0.002; 0.040; 0.078; 0.115; 0.150; 0.182; 0.212; 0.240; 0.265; 0.288; 0.308; 0.325; 0.340];

assay_bank = [assay_control_S1, assay_control_S2, assay_inhibitor_S2];
n_assays = size(assay_bank, 2);

v0_calculated = zeros(n_assays, 1);
window_points = 5; % Minimum points to evaluate the moving linearity window

figure('Name', 'Module I: Linear-Phase Diagnosis (v0)', 'Color', 'w', 'Position', [100, 500, 550, 400]);

for j = 1:n_assays
    current_signal = assay_bank(:, j);
    best_R2 = 0;
    best_slope = 0;
    idx_start = 1; idx_end = window_points;

    % Moving-window algorithm
    for i = 1 : (length(common_time) - window_points + 1)
        window = i : (i + window_points - 1);
        t_sub = common_time(window);
        s_sub = current_signal(window);
        p = polyfit(t_sub, s_sub, 1);
        R2 = 1 - (sum((s_sub - polyval(p, t_sub)).^2) / sum((s_sub - mean(s_sub)).^2));
        if R2 > best_R2 && i <= (length(common_time)/2)
            best_R2 = R2;
            best_slope = p(1);
            idx_start = i; idx_end = i + window_points - 1;
        end
    end
    v0_calculated(j) = best_slope; % Store the calculated v0

    % Plot linear-phase quality control
    plot(common_time, current_signal, 'o-', 'LineWidth', 1); hold on;
    t_lin = common_time(idx_start:idx_end);
    plot(t_lin, polyval(polyfit(t_lin, current_signal(idx_start:idx_end), 1), t_lin), 'k-', 'LineWidth', 2);
end
xlabel('Time (s)'); ylabel('Signal / Absorbance');
title('Module I: Automated v0 Extraction'); grid on;
fprintf('>> MODULE I COMPLETED: v0 values extracted successfully.\n\n');

%% =========================================================================
% MODULE II: MACROSCOPIC, MULTI-MODEL ANALYSIS AND INHIBITION DIAGNOSIS
% =========================================================================
disp('>> MODULE II: Running macroscopic modeling and AIC calculation...');

% MAIN ENZYMATIC WORK MATRIX
% You can feed this matrix with the results of Module I.
% Format: [ Substrate_Concentration [S] , Inhibitor_Concentration [I] , Initial_Velocity v ]
% (Integrated example of a complete catalase-style kinetic assay at high
% substrate concentration)
macro_data = [
%   [S]    [I]     v
    5.0,   0.0,   45.0;
    10.0,  0.0,   78.0;
    20.0,  0.0,  120.0;
    40.0,  0.0,  145.0; % Typical curve decaying by substrate inhibition
    60.0,  0.0,  138.0;
    80.0,  0.0,  122.0;
    100.0, 0.0,  105.0;
];

S = macro_data(:, 1);
I = macro_data(:, 2);
v = macro_data(:, 3);
N = length(v);

% Intelligently detect whether an external inhibitor is present in the assay
has_inhibitor = any(I > 0);

% GENERAL BANK OF ANONYMOUS MATHEMATICAL MODELS
mod_mm      = @(beta, x) (beta(1) * x) ./ (beta(2) + x);                                          % Michaelis-Menten
mod_hill    = @(beta, x) (beta(1) * x.^beta(3)) ./ (beta(2)^beta(3) + x.^beta(3));                % Hill (empirical allostery)
mod_haldane = @(beta, x) (beta(1) * x) ./ (beta(2) + x + (x.^2 / beta(3)));                       % Haldane (catalase / substrate inhibition)
mod_act_sust= @(beta, x) (beta(1) * x .* (1 + (beta(4)*x)/beta(3))) ./ (beta(2) + x .* (1 + x/beta(3))); % Substrate activation
mod_mwc     = @(beta, x) beta(1) * ( (x./beta(2)).*(1 + x./beta(2)) + beta(3)*beta(4)*(x./beta(2)).*(1 + beta(4)*x./beta(2)) ) ./ ...
                  ( (1 + x./beta(2)).^2 + beta(3)*(1 + beta(4)*x./beta(2)).^2 );                  % MWC (mechanistic)
mod_comp    = @(beta, x) (beta(1) * x(:,1)) ./ (x(:,1) + beta(2) * (1 + x(:,2)/beta(3)));         % Competitive inhibition
mod_nocomp  = @(beta, x) (beta(1) * x(:,1)) ./ ((beta(2) + x(:,1)) .* (1 + x(:,2)/beta(3)));      % Pure non-competitive inhibition
mod_incomp  = @(beta, x) (beta(1) * x(:,1)) ./ (beta(2) + x(:,1) .* (1 + x(:,2)/beta(3)));        % Uncompetitive inhibition
mod_mixed   = @(beta, x) (beta(1) * x(:,1)) ./ (beta(2)*(1 + x(:,2)/beta(3)) + x(:,1).*(1 + x(:,2)/(beta(4)*beta(3)))); % Mixed inhibition

% Conditional execution according to the type of input data
if ~has_inhibitor
    % ---------------------------------------------------------------------
    % SCENARIO WITHOUT INHIBITOR: NON-MICHAELIAN COMPETITION
    % ---------------------------------------------------------------------
    disp('>> Active analysis: substrate-dependent scenario (no external inhibitor).');

    % Nonlinear fits
    p0_mm = [max(v)*2, median(S)];
    [p_mm, res_mm] = lsqcurvefit(mod_mm, p0_mm, S, v, [0,0], [Inf,Inf], options);
    AIC_mm = N * log(res_mm/N) + 2 * length(p0_mm);

    p0_hill = [max(v), median(S), 1.0];
    [p_hill, res_hill] = lsqcurvefit(mod_hill, p0_hill, S, v, [0,0,0], [Inf,Inf,10], options);
    AIC_hill = N * log(res_hill/N) + 2 * length(p0_hill);

    p0_hald = [max(v)*2, median(S), max(S)];
    [p_hald, res_hald] = lsqcurvefit(mod_haldane, p0_hald, S, v, [0,0,0], [Inf,Inf,Inf], options);
    AIC_hald = N * log(res_hald/N) + 2 * length(p0_hald);

    p0_mwc = [max(v), median(S), 10, 0.1];
    [p_mwc, res_mwc] = lsqcurvefit(mod_mwc, p0_mwc, S, v, [0,0,0,0], [Inf,Inf,Inf,1], options);
    AIC_mwc = N * log(res_mwc/N) + 2 * length(p0_mwc);

    p0_act = [max(v), median(S), median(S), 2.0];
    [p_act, res_act] = lsqcurvefit(mod_act_sust, p0_act, S, v, [0,0,0,1], [Inf,Inf,Inf,Inf], options);
    AIC_act = N * log(res_act/N) + 2 * length(p0_act);

    Models      = {'Michaelis-Menten'; 'Hill (Empirical Allosteric)'; 'Haldane (Catalase/Substrate Inhibition)'; 'MWC (Mechanistic Allostery)'; 'Substrate Activation'};
    Residuals_SS = [res_mm; res_hill; res_hald; res_mwc; res_act];
    AIC_values   = [AIC_mm; AIC_hill; AIC_hald; AIC_mwc; AIC_act];
else
    % ---------------------------------------------------------------------
    % SCENARIO WITH INHIBITOR: COMPETITION OF INHIBITION TYPES
    % ---------------------------------------------------------------------
    disp('>> Active analysis: scenario with external inhibitor.');

    p0_3par = [max(v)*1.2, median(S), max(I)];
    p0_4par = [max(v)*1.2, median(S), max(I), 1.0];

    [p_comp, res_comp] = lsqcurvefit(mod_comp, p0_3par, [S, I], v, [0,0,0], [Inf,Inf,Inf], options);
    AIC_comp = N * log(res_comp/N) + 2 * length(p0_3par);

    [p_nocomp, res_nocomp] = lsqcurvefit(mod_nocomp, p0_3par, [S, I], v, [0,0,0], [Inf,Inf,Inf], options);
    AIC_nocomp = N * log(res_nocomp/N) + 2 * length(p0_3par);

    [p_incomp, res_incomp] = lsqcurvefit(mod_incomp, p0_3par, [S, I], v, [0,0,0], [Inf,Inf,Inf], options);
    AIC_incomp = N * log(res_incomp/N) + 2 * length(p0_3par);

    [p_mixed, res_mixed] = lsqcurvefit(mod_mixed, p0_4par, [S, I], v, [0,0,0,0], [Inf,Inf,Inf,Inf], options);
    AIC_mixed = N * log(res_mixed/N) + 2 * length(p0_4par);

    Models       = {'Competitive Inhibition'; 'Pure Non-Competitive Inhibition'; 'Uncompetitive Inhibition'; 'Mixed Inhibition'};
    Residuals_SS = [res_comp; res_nocomp; res_incomp; res_mixed];
    AIC_values   = [AIC_comp; AIC_nocomp; AIC_incomp; AIC_mixed];
end

% =========================================================================
% FINAL REPORT OF CONCLUSIONS AND KINETIC PARAMETERS
% =========================================================================
ResultsTable = table(Models, Residuals_SS, AIC_values);
disp(' '); disp('---------------------------------------------------------');
disp(' COMPARATIVE MATHEMATICAL MODELING TABLE ');
disp('---------------------------------------------------------');
disp(ResultsTable);

[~, winner_idx] = min(AIC_values);
fprintf('*** SCIENTIFIC VERDICT: The optimal model is [%s] ***\n\n', Models{winner_idx});
fprintf('--- KINETIC PARAMETERS EXTRACTED FROM THE EXPERIMENT ---\n');

if ~has_inhibitor
    switch winner_idx
        case 1 % MM
            fprintf('Vmax = %.2f\nKm = %.2f\n', p_mm(1), p_mm(2));
        case 2 % Hill
            fprintf('Vmax = %.2f\nK_half = %.2f\nHill coefficient (n) = %.2f\n', p_hill(1), p_hill(2), p_hill(3));
        case 3 % Haldane
            fprintf('Theoretical Vmax (free enzyme) = %.2f\nApparent Km = %.2f\nKsi (substrate self-inhibition constant) = %.2f\n', p_hald(1), p_hald(2), p_hald(3));
            fprintf('>> Calculated optimal substrate concentration: %.2f\n', sqrt(p_hald(2)*p_hald(3)));
        case 4 % MWC
            fprintf('Vmax = %.2f\nKR = %.2f\nL (T/R equilibrium) = %.2f\nc (KR/KT ratio) = %.2f\n', p_mwc(1), p_mwc(2), p_mwc(3), p_mwc(4));
        case 5 % Substrate activation
            fprintf('Basal Vmax = %.2f\nKm = %.2f\nKA (activation) = %.2f\nBeta acceleration factor = %.2f\n', p_act(1), p_act(2), p_act(3), p_act(4));
    end
else
    switch winner_idx
        case 1 % Competitive
            fprintf('Vmax = %.2f\nKm = %.2f\nKi (inhibitor affinity) = %.2f\n', p_comp(1), p_comp(2), p_comp(3));
        case 2 % Non-competitive
            fprintf('Vmax = %.2f\nKm = %.2f\nKi (inhibitor affinity) = %.2f\n', p_nocomp(1), p_nocomp(2), p_nocomp(3));
        case 3 % Uncompetitive
            fprintf('Vmax = %.2f\nKm = %.2f\nKi (affinity for the ES complex) = %.2f\n', p_incomp(1), p_incomp(2), p_incomp(3));
        case 4 % Mixed
            fprintf('Vmax = %.2f\nKm = %.2f\nKi (free site) = %.2f\nAlpha*Ki (occupied site) = %.2f\nAlpha factor = %.2f\n', p_mixed(1), p_mixed(2), p_mixed(3), p_mixed(4)*p_mixed(3), p_mixed(4));
    end
end

% =========================================================================
% PLOTTING OF FITS IN THE MACROSCOPIC MODELING
% =========================================================================
S_plot = linspace(0, max(S)*1.2, 300);
figure('Name', 'Module II: Macroscopic Fitting Curves', 'Color', 'w', 'Position', [700, 500, 600, 400]);
plot(S, v, 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Real v0 data'); hold on;

if ~has_inhibitor
    plot(S_plot, mod_mm(p_mm, S_plot), 'g--', 'DisplayName', 'Michaelis-Menten');
    plot(S_plot, mod_haldane(p_hald, S_plot), 'b-', 'LineWidth', 1.8, 'DisplayName', 'Haldane (Winner)');
    plot(S_plot, mod_mwc(p_mwc, S_plot), 'm-.', 'DisplayName', 'MWC Model');
    plot(S_plot, mod_act_sust(p_act, S_plot), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Substrate Activation');
    title('Non-Michaelian Kinetics Fitting');
else
    I_max = max(I);
    if winner_idx == 1
        v_ctrl = mod_comp(p_comp,[S_plot',zeros(300,1)]); v_inh = mod_comp(p_comp,[S_plot',ones(300,1)*I_max]);
    elseif winner_idx == 2
        v_ctrl = mod_nocomp(p_nocomp,[S_plot',zeros(300,1)]); v_inh = mod_nocomp(p_nocomp,[S_plot',ones(300,1)*I_max]);
    elseif winner_idx == 3
        v_ctrl = mod_incomp(p_incomp,[S_plot',zeros(300,1)]); v_inh = mod_incomp(p_incomp,[S_plot',ones(300,1)*I_max]);
    else
        v_ctrl = mod_mixed(p_mixed,[S_plot',zeros(300,1)]); v_inh = mod_mixed(p_mixed,[S_plot',ones(300,1)*I_max]);
    end
    plot(S_plot, v_ctrl, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Control fit');
    plot(S_plot, v_inh, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Fit with inhibitor');
    title('Enzyme Inhibition Profile');
end
xlabel('Substrate Concentration [S]'); ylabel('Initial Velocity (v0)');
legend('Location', 'Best'); grid on;
