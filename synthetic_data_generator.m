%% synthetic_data_generator.m
%% SYNTHETIC DATASET GENERATOR FOR ENZYME KINETICS BENCHMARKING
% Generates synthetic initial-velocity datasets with known ground truth:
% the user selects any model implemented in the SSICE catalog, assigns
% "true" parameters, and adds adjustable Gaussian noise. Because the
% generating model and parameters are known, the datasets are used to
% benchmark parameter recovery and AICc-based model discrimination under
% controlled conditions (see Section 2.8 of the SSICE manuscript).
%
% The output matrix [S, I, v] can be pasted directly into the macro_data
% matrix of ssice_integrated_kinetics.m (Module II), or saved as a
% two/three-column file loadable from michaelis_menten_analysis.m.
%
% Requirements: MATLAB (base installation; no toolboxes required)
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

clear; clc; close all;
fprintf('=========================================================================\n');
fprintf(' SYNTHETIC DATASET GENERATOR - KNOWN GROUND TRUTH FOR BENCHMARKING\n');
fprintf('=========================================================================\n\n');

%% =========================================================================
% 1. CONFIGURATION (edit this block)
% =========================================================================
% Model selection:
%   Substrate-dependent models (no inhibitor):
%     'michaelis-menten' | 'hill' | 'haldane' | 'mwc' | 'substrate-activation'
%   Inhibition models (require I_levels below):
%     'competitive' | 'noncompetitive' | 'uncompetitive' | 'mixed'
model_choice = 'haldane';

% True parameters (defaults within the realistic ranges of Section 2.8:
% Vmax 0.01-1 mM/s, Km 0.1-10 mM). Only the parameters of the selected
% model are used:
%   michaelis-menten     : [Vmax, Km]
%   hill                 : [Vmax, K_half, n]
%   haldane              : [Vmax, Km, Ksi]          (substrate inhibition)
%   mwc                  : [Vmax, KR, L, c]
%   substrate-activation : [Vmax, Km, KA, beta]
%   competitive          : [Vmax, Km, Ki]
%   noncompetitive       : [Vmax, Km, Ki]
%   uncompetitive        : [Vmax, Km, Ki]
%   mixed                : [Vmax, Km, Ki, alpha]
true_params = [0.5, 5.0, 250.0];          % Haldane example

% Experimental design
n_points  = 10;                            % number of substrate levels
S_min     = 0.5;                           % lowest [S] (mM)
S_max     = 100.0;                         % highest [S] (mM)
I_levels  = [0, 5, 10, 20, 50];            % inhibitor levels (only for
                                           % inhibition models)
% Gaussian noise: coefficient of variation (CV) as a fraction of the rate
% value. Manuscript benchmark: 0.05 (sigma = 5% of v); adjustable 0.01-0.20.
noise_cv  = 0.05;

% Reproducibility and output
random_seed = 42;                          % rng seed for reproducibility
save_file   = 'synthetic_dataset.csv';     % '' to skip file export
make_figure = true;

%% =========================================================================
% 2. MODEL BANK (identical functional forms as ssice_integrated_kinetics.m)
% =========================================================================
models = struct();
models.michaelis_menten = @(b, S, I) (b(1) .* S) ./ (b(2) + S);
models.hill             = @(b, S, I) (b(1) .* S.^b(3)) ./ (b(2)^b(3) + S.^b(3));
models.haldane          = @(b, S, I) (b(1) .* S) ./ (b(2) + S + (S.^2 ./ b(3)));
models.mwc              = @(b, S, I) b(1) .* ( (S./b(2)).*(1 + S./b(2)) + b(3).*b(4).*(S./b(2)).*(1 + b(4).*S./b(2)) ) ./ ...
                              ( (1 + S./b(2)).^2 + b(3).*(1 + b(4).*S./b(2)).^2 );
models.substrate_activation = @(b, S, I) (b(1) .* S .* (1 + (b(4).*S)/b(3))) ./ (b(2) + S .* (1 + S./b(3)));
models.competitive      = @(b, S, I) (b(1) .* S) ./ (S + b(2) .* (1 + I./b(3)));
models.noncompetitive   = @(b, S, I) (b(1) .* S) ./ ((b(2) + S) .* (1 + I./b(3)));
models.uncompetitive    = @(b, S, I) (b(1) .* S) ./ (b(2) + S .* (1 + I./b(3)));
models.mixed            = @(b, S, I) (b(1) .* S) ./ (b(2).*(1 + I./b(3)) + S.*(1 + I./(b(4).*b(3))));

model_field = strrep(model_choice, '-', '_');
if ~isfield(models, model_field)
    error('Unknown model "%s". See the model list in the configuration block.', model_choice);
end
model_fun = models.(model_field);

% Inhibition models require inhibitor levels; substrate-dependent models
% are generated with a single I = 0 column for matrix compatibility.
inhibition_models = {'competitive', 'noncompetitive', 'uncompetitive', 'mixed'};
is_inhibition = ismember(model_choice, inhibition_models);
if ~is_inhibition
    I_levels = 0;
end

%% =========================================================================
% 3. DATASET GENERATION
% =========================================================================
rng(random_seed);

% Substrate grid (log-spaced to sample low and high [S] evenly)
S_grid = logspace(log10(S_min), log10(S_max), n_points)';

n_S = length(S_grid);
n_I = length(I_levels);
data = zeros(n_S * n_I, 3);   % columns: [S, I, v_observed]

fprintf('>> Generating dataset with model: %s\n', model_choice);
fprintf('   True parameters: %s\n', mat2str(true_params, 4));
fprintf('   Gaussian noise: CV = %.0f %% of the rate value\n\n', noise_cv * 100);

row = 1;
for I = I_levels
    % Exact (noise-free) velocities from the ground-truth model
    v_true = model_fun(true_params, S_grid, I * ones(n_S, 1));

    % Proportional Gaussian noise: sigma = CV * v
    noise  = noise_cv .* v_true .* randn(n_S, 1);
    v_obs  = v_true + noise;
    v_obs  = max(v_obs, 0);   % physical constraint: no negative rates

    data(row:row+n_S-1, :) = [S_grid, I * ones(n_S, 1), v_obs];
    row = row + n_S;
end

%% =========================================================================
% 4. GROUND-TRUTH REPORT (for parameter-recovery benchmarking)
% =========================================================================
fprintf('--- GROUND TRUTH (keep this table for the recovery benchmark) ---\n');
param_names = param_name_list(model_choice);
for k = 1:numel(param_names)
    fprintf('   %-8s = %.6g\n', param_names{k}, true_params(k));
end
fprintf('   Model    = %s\n', model_choice);
fprintf('   Noise CV = %.0f %%\n', noise_cv * 100);
fprintf('   Points   = %d substrate levels x %d inhibitor levels = %d rows\n\n', ...
    n_S, n_I, size(data, 1));
fprintf(['>> Benchmark protocol: (1) feed this dataset into Module II of\n', ...
         '   ssice_integrated_kinetics.m without disclosing the model;\n', ...
         '   (2) check that AICc ranks "%s" first (model discrimination);\n', ...
         '   (3) compare the fitted parameters against the ground truth\n', ...
         '   above (parameter recovery).\n\n'], model_choice);

%% =========================================================================
% 5. FILE EXPORT (three columns: [S]  [I]  v)
% =========================================================================
if ~isempty(save_file)
    writematrix(data, save_file);
    fprintf('>> Dataset saved to: %s\n\n', save_file);
end

%% =========================================================================
% 6. DIAGNOSTIC FIGURE (true curve vs. noisy observations)
% =========================================================================
if make_figure
    figure('Name', 'Synthetic Dataset - Ground Truth vs. Observations', ...
        'Color', 'w', 'Position', [150, 150, 650, 420]);
    S_fine = logspace(log10(S_min), log10(S_max * 1.1), 300)';
    cmap = parula(n_I);
    hold on;
    for idx = 1:n_I
        I = I_levels(idx);
        v_curve = model_fun(true_params, S_fine, I * ones(300, 1));
        rows = data(:, 2) == I;
        plot(S_fine, v_curve, '-', 'Color', cmap(idx, :), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('True curve, [I] = %g', I));
        plot(data(rows, 1), data(rows, 3), 'o', 'Color', cmap(idx, :), ...
            'MarkerFaceColor', cmap(idx, :), 'MarkerSize', 5, ...
            'HandleVisibility', 'off');
    end
    set(gca, 'XScale', 'log');
    xlabel('Substrate concentration [S] (mM)', 'FontSize', 11);
    ylabel('Initial velocity v', 'FontSize', 11);
    title(sprintf('Synthetic dataset: %s (CV = %.0f%%)', ...
        strrep(model_choice, '_', '-'), noise_cv * 100), 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 9);
    grid on;
end

%% =========================================================================
% AUXILIARY FUNCTION
% =========================================================================
function names = param_name_list(model_choice)
    switch model_choice
        case 'michaelis-menten'
            names = {'Vmax', 'Km'};
        case 'hill'
            names = {'Vmax', 'K_half', 'n'};
        case 'haldane'
            names = {'Vmax', 'Km', 'Ksi'};
        case 'mwc'
            names = {'Vmax', 'KR', 'L', 'c'};
        case 'substrate-activation'
            names = {'Vmax', 'Km', 'KA', 'beta'};
        case {'competitive', 'noncompetitive', 'uncompetitive'}
            names = {'Vmax', 'Km', 'Ki'};
        case 'mixed'
            names = {'Vmax', 'Km', 'Ki', 'alpha'};
        otherwise
            names = {};
    end
end
