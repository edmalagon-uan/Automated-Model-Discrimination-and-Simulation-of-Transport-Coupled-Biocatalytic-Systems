%% michaelis_menten_analysis.m
% MAIN PROGRAM FOR ENZYME KINETICS ANALYSIS (MICHAELIS-MENTEN)
%
% Interactive tool to estimate Vmax and Km from initial-velocity data
% using (i) direct nonlinear Michaelis-Menten fitting and (ii) the three
% classical linear transformations (Lineweaver-Burk, Eadie-Hofstee and
% Hanes-Woolf). Produces the multi-panel diagnostic figure used as
% Figure 1 of the SSICE manuscript.
%
% Requirements: MATLAB R2016b or later (local functions in scripts)
%               Optimization Toolbox (lsqcurvefit)
%
% Usage: run the script and choose a data source from the main menu:
%   1. Use example data
%   2. Enter data manually
%   3. Load data from file (two columns: [S] and v)
%   4. Simulate data with parameter variation
%   5. Exit
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

clear all; close all; clc;
fprintf('============================================\n');
fprintf('ENZYME KINETICS ANALYSIS - MICHAELIS-MENTEN\n');
fprintf('============================================\n\n');

% MAIN MENU
option = 0;
while option ~= 5
    fprintf('\n--- MAIN MENU ---\n');
    fprintf('1. Use example data\n');
    fprintf('2. Enter data manually\n');
    fprintf('3. Load data from file\n');
    fprintf('4. Simulate data with parameter variation\n');
    fprintf('5. Exit\n');
    option = input('Select an option (1-5): ');
    switch option
        case 1
            example_data();
        case 2
            enter_data();
        case 3
            load_file();
        case 4
            simulate_parameters();
        case 5
            fprintf('Goodbye!\n');
        otherwise
            fprintf('Invalid option. Please try again.\n');
    end
end

%% FUNCTION FOR EXAMPLE DATA
function example_data()
    fprintf('\n--- USING EXAMPLE DATA ---\n');
    % Example data (substrate concentration vs. velocity)
    S = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]; % mM
    v = [2.1, 3.8, 6.5, 8.2, 9.1, 9.8, 10.2, 10.4]; % umol/min
    % Analyze data
    analyze_kinetics(S, v);
end

%% FUNCTION FOR MANUAL DATA ENTRY
function enter_data()
    fprintf('\n--- MANUAL DATA ENTRY ---\n');
    n = input('Number of experimental points: ');
    S = zeros(1, n);
    v = zeros(1, n);
    fprintf('\nEnter the [S] (mM) and v (umol/min) pairs:\n');
    for i = 1:n
        fprintf('Point %d:\n', i);
        S(i) = input('  [S] (mM): ');
        v(i) = input('  v (umol/min): ');
    end
    % Analyze data
    analyze_kinetics(S, v);
end

%% FUNCTION TO LOAD DATA FROM FILE
function load_file()
    fprintf('\n--- LOAD DATA FROM FILE ---\n');
    fprintf('The file must contain two columns: [S] and v\n');
    file_name = input('File name (with extension): ', 's');
    try
        data = load(file_name);
        S = data(:, 1);
        v = data(:, 2);
        fprintf('Data loaded successfully:\n');
        fprintf('  %d experimental points\n', length(S));
        % Analyze data
        analyze_kinetics(S, v);
    catch
        fprintf('Error loading the file. Check the name and format.\n');
    end
end

%% FUNCTION TO SIMULATE WITH PARAMETER VARIATION
function simulate_parameters()
    fprintf('\n--- SIMULATION WITH PARAMETER VARIATION ---\n');
    % Simulation parameters
    Vmax_true = input('True value of Vmax (umol/min): ');
    Km_true   = input('True value of Km (mM): ');
    noise     = input('Noise level (0 = no noise, 1 = maximum): ');
    % Generate simulated data
    S_sim = logspace(-1, 2, 15); % Concentrations on a logarithmic scale
    v_sim = (Vmax_true * S_sim) ./ (Km_true + S_sim);
    % Add noise if requested
    if noise > 0
        random_noise = noise * 0.1 * Vmax_true * randn(size(v_sim));
        v_sim = v_sim + random_noise;
        v_sim = max(v_sim, 0); % Do not allow negative velocities
    end
    % Analyze the simulated data
    fprintf('\nAnalyzing simulated data...\n');
    [Vmax_est, Km_est, ~] = analyze_kinetics(S_sim, v_sim);
    % Compare with true values
    fprintf('\n--- PARAMETER COMPARISON ---\n');
    fprintf('Parameter      True     Estimated   Error (%%)\n');
    fprintf('Vmax:     %10.2f %10.2f %10.2f\n', ...
        Vmax_true, Vmax_est, abs((Vmax_est - Vmax_true)/Vmax_true*100));
    fprintf('Km:       %10.2f %10.2f %10.2f\n', ...
        Km_true, Km_est, abs((Km_est - Km_true)/Km_true*100));
end

%% MAIN KINETIC ANALYSIS FUNCTION
function [Vmax, Km, stats] = analyze_kinetics(S, v)
    fprintf('\n--- ENZYME KINETICS ANALYSIS ---\n');

    % 1. Direct Michaelis-Menten fit (nonlinear)
    fprintf('\n1. Nonlinear fit (direct Michaelis-Menten):\n');
    % Initial guess
    Vmax_guess = max(v);
    Km_guess = S(find(v >= Vmax_guess/2, 1));
    % Michaelis-Menten model
    mm_model = @(params, S) (params(1) * S) ./ (params(2) + S);
    % Nonlinear fit
    params0 = [Vmax_guess, Km_guess];
    params_opt = lsqcurvefit(mm_model, params0, S, v);
    Vmax_nl = params_opt(1);
    Km_nl = params_opt(2);
    fprintf('   Vmax = %.4f umol/min\n', Vmax_nl);
    fprintf('   Km   = %.4f mM\n', Km_nl);

    % 2. Lineweaver-Burk transformation (double reciprocal)
    fprintf('\n2. Lineweaver-Burk transformation:\n');
    % Compute reciprocals
    x_lb = 1./S;
    y_lb = 1./v;
    % Linear fit
    p_lb = polyfit(x_lb, y_lb, 1);
    Vmax_lb = 1/p_lb(2);
    Km_lb = p_lb(1)/p_lb(2);
    fprintf('   Vmax = %.4f umol/min\n', Vmax_lb);
    fprintf('   Km   = %.4f mM\n', Km_lb);

    % 3. Eadie-Hofstee transformation
    fprintf('\n3. Eadie-Hofstee transformation:\n');
    x_eh = v./S;
    y_eh = v;
    p_eh = polyfit(x_eh, y_eh, 1);
    Vmax_eh = p_eh(2);
    Km_eh = -p_eh(1);
    fprintf('   Vmax = %.4f umol/min\n', Vmax_eh);
    fprintf('   Km   = %.4f mM\n', Km_eh);

    % 4. Hanes-Woolf transformation
    fprintf('\n4. Hanes-Woolf transformation:\n');
    x_hw = S;
    y_hw = S./v;
    p_hw = polyfit(x_hw, y_hw, 1);
    Vmax_hw = 1/p_hw(1);
    Km_hw = p_hw(2)/p_hw(1);
    fprintf('   Vmax = %.4f umol/min\n', Vmax_hw);
    fprintf('   Km   = %.4f mM\n', Km_hw);

    % Final parameters (nonlinear method)
    Vmax = Vmax_nl;
    Km = Km_nl;

    % Compute statistics
    v_pred = (Vmax * S) ./ (Km + S);
    residuals = v - v_pred;
    SSE = sum(residuals.^2);
    SST = sum((v - mean(v)).^2);
    R2 = 1 - SSE/SST;
    RMSE = sqrt(mean(residuals.^2));
    stats = struct('SSE', SSE, 'R2', R2, 'RMSE', RMSE);
    fprintf('\n--- FIT STATISTICS ---\n');
    fprintf('R^2  = %.4f\n', R2);
    fprintf('RMSE = %.4f umol/min\n', RMSE);
    fprintf('SSE  = %.4f\n', SSE);

    % PLOTS
    create_plots(S, v, v_pred, Vmax, Km, x_lb, y_lb, p_lb);
end

%% FUNCTION TO CREATE PLOTS
function create_plots(S, v, v_pred, Vmax, Km, x_lb, y_lb, p_lb)
    % Figure 1: Enzyme kinetics analysis
    figure('Name', 'Enzyme Kinetics Analysis', 'Position', [100, 100, 1200, 800]);

    % Subplot 1: Michaelis-Menten
    subplot(2, 3, 1);
    S_fine = linspace(0, max(S)*1.2, 100);
    v_fine = (Vmax * S_fine) ./ (Km + S_fine);
    plot(S, v, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    plot(S_fine, v_fine, 'r-', 'LineWidth', 2);
    plot([0, Km], [Vmax/2, Vmax/2], 'g--', 'LineWidth', 1);
    plot([Km, Km], [0, Vmax/2], 'g--', 'LineWidth', 1);
    xlabel('[S] (mM)', 'FontSize', 12);
    ylabel('v (umol/min)', 'FontSize', 12);
    title('Michaelis-Menten', 'FontSize', 14);
    legend('Data', 'Fit', 'V_{max}/2', 'K_m', 'Location', 'southeast');
    grid on;

    % Subplot 2: Lineweaver-Burk
    subplot(2, 3, 2);
    x_fit = linspace(min(x_lb)*0.9, max(x_lb)*1.1, 50);
    y_fit = polyval(p_lb, x_fit);
    plot(x_lb, y_lb, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
    xlabel('1/[S] (mM^{-1})', 'FontSize', 12);
    ylabel('1/v (min/umol)', 'FontSize', 12);
    title('Lineweaver-Burk', 'FontSize', 14);
    legend('Data', 'Fit', 'Location', 'best');
    grid on;

    % Subplot 3: Eadie-Hofstee
    subplot(2, 3, 3);
    v_s = v./S;
    plot(v_s, v, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    x_eh_fit = linspace(min(v_s), max(v_s), 50);
    y_eh_fit = Vmax - Km * x_eh_fit;
    plot(x_eh_fit, y_eh_fit, 'r-', 'LineWidth', 2);
    xlabel('v/[S] (min^{-1})', 'FontSize', 12);
    ylabel('v (umol/min)', 'FontSize', 12);
    title('Eadie-Hofstee', 'FontSize', 14);
    legend('Data', 'Fit', 'Location', 'best');
    grid on;

    % Subplot 4: Hanes-Woolf
    subplot(2, 3, 4);
    S_v = S./v;
    plot(S, S_v, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    p_hw = polyfit(S, S_v, 1);
    x_hw_fit = linspace(0, max(S), 50);
    y_hw_fit = polyval(p_hw, x_hw_fit);
    plot(x_hw_fit, y_hw_fit, 'r-', 'LineWidth', 2);
    xlabel('[S] (mM)', 'FontSize', 12);
    ylabel('[S]/v (min·mM/umol)', 'FontSize', 12);
    title('Hanes-Woolf', 'FontSize', 14);
    legend('Data', 'Fit', 'Location', 'best');
    grid on;

    % Subplot 5: Residuals
    subplot(2, 3, 5);
    residuals = v - v_pred;
    plot(S, residuals, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    hold on;
    plot([0, max(S)], [0, 0], 'r--', 'LineWidth', 1);
    xlabel('[S] (mM)', 'FontSize', 12);
    ylabel('Residuals (umol/min)', 'FontSize', 12);
    title('Residual Analysis', 'FontSize', 14);
    grid on;

    % Subplot 6: Information
    subplot(2, 3, 6);
    axis off;
    text(0.1, 0.9, 'ANALYSIS RESULTS', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('V_{max} = %.3f umol/min', Vmax), 'FontSize', 12);
    text(0.1, 0.6, sprintf('K_m = %.3f mM', Km), 'FontSize', 12);
    text(0.1, 0.4, 'Equation:', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.25, sprintf('v = %.3f [S] / (%.3f + [S])', Vmax, Km), 'FontSize', 11);
    text(0.1, 0.1, 'Units: [S] in mM, v in umol/min', 'FontSize', 10, 'FontStyle', 'italic');

    % Additional figure: comparison of methods
    figure('Name', 'Comparison of Methods', 'Position', [150, 150, 800, 600]);
    % Compute predictions of each method
    S_fine = linspace(0, max(S)*1.2, 100);
    v_mm = (Vmax * S_fine) ./ (Km + S_fine);
    plot(S, v, 'ko', 'MarkerSize', 10, 'LineWidth', 2);
    hold on;
    plot(S_fine, v_mm, 'r-', 'LineWidth', 3);
    xlabel('[S] (mM)', 'FontSize', 14);
    ylabel('v (umol/min)', 'FontSize', 14);
    title('Comparison of Fitting Methods', 'FontSize', 16);
    legend('Experimental data', 'Michaelis-Menten (nonlinear)', 'Location', 'southeast');
    grid on;
    % Show parameters on the plot
    text(max(S)*0.6, max(v)*0.9, sprintf('V_{max} = %.3f', Vmax), 'FontSize', 12);
    text(max(S)*0.6, max(v)*0.8, sprintf('K_m = %.3f', Km), 'FontSize', 12);
    text(max(S)*0.6, max(v)*0.7, sprintf('R^2 = %.3f', 1 - sum((v - v_pred).^2)/sum((v - mean(v)).^2)), 'FontSize', 12);
end
