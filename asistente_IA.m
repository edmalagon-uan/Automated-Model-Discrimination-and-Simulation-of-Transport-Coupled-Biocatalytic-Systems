%% asistente_IA.m
%% OPTIONAL LLM-BASED ASSISTANT FOR THE SSICE FRAMEWORK
% Builds structured prompts from the numerical results of the kinetic
% pipeline (candidate models, AICc values, parameter estimates with
% confidence intervals, residual statistics, R^2) and returns
% natural-language guidance templates for four stages:
%   (i)   data inspection
%   (ii)  model pre-selection from curve shape
%   (iii) post-fit diagnosis (AICc differences, wide confidence
%         intervals, systematic residuals)
%   (iv)  simulation guidance
%
% SAFEGUARDS (see Section 2.7 of the SSICE manuscript):
%   1. A hard-coded heuristic layer produces a rule-based pre-diagnosis
%      against which the LLM recommendations must be cross-checked.
%   2. All suggestions are presented as hypotheses to be verified by
%      rerunning the numerical analysis.
%   3. Every output is prefixed by a disclaimer: the assistant does not
%      replace expert judgment or the numerical results.
%
% The assistant NEVER modifies scripts, parameters, or results; only
% aggregated statistics leave the pipeline (in the prompt text).
%
% USAGE:
%   results = struct with fields (see build_example_results below):
%       models, AICc, params, param_names, ci, R2, residuals, S, v
%   stage = 1..4
%   asistente_IA(results, stage)
%
% The script prints (and optionally saves) the prompt to paste into any
% LLM interface. After the LLM session, run crosscheck_response() with
% the model's recommendation to document the safeguard agreement.
% Record the model name, version, and session date for the Supporting
% Information (Section S1).
%
% Requirements: MATLAB R2016b or later (local functions in scripts);
%               no external toolboxes and no API keys required.
%
% Part of the SSICE framework. If you use this code, please cite the
% associated manuscript (see README.md).

function asistente_IA(results, stage, save_prompt)
    if nargin < 3, save_prompt = ''; end

    fprintf('=========================================================================\n');
    fprintf(' SSICE LLM ASSISTANT - PROMPT BUILDER AND SAFEGUARD LAYER\n');
    fprintf('=========================================================================\n\n');
    fprintf('DISCLAIMER: The guidance produced with this assistant consists of\n');
    fprintf('hypotheses to be verified by rerunning the numerical analysis. It does\n');
    fprintf('not replace expert judgment or the numerical results.\n\n');

    % 1. Structured summary of the numerical results
    summary = summarize_results(results);
    fprintf('--- STRUCTURED SUMMARY SENT TO THE ASSISTANT ---\n%s\n', summary);

    % 2. Heuristic safeguard layer (rule-based pre-diagnosis)
    flags = heuristic_rules(results);
    fprintf('--- HEURISTIC PRE-DIAGNOSIS (safeguard layer) ---\n');
    if isempty(flags)
        fprintf('   No rule-based warnings triggered.\n\n');
    else
        for k = 1:numel(flags)
            fprintf('   [%d] %s\n', k, flags{k});
        end
        fprintf('\n');
    end

    % 3. Build the stage-specific prompt
    prompt = build_prompt(summary, flags, stage);
    fprintf('--- PROMPT READY (paste into your LLM interface) ---\n%s\n', prompt);

    if ~isempty(save_prompt)
        fid = fopen(save_prompt, 'w');
        fprintf(fid, '%s\n', prompt);
        fclose(fid);
        fprintf('>> Prompt saved to: %s\n', save_prompt);
    end
    fprintf(['>> After the LLM session: record the model name, version,\n' ...
             '   and date (needed for the Supporting Information, S1), and\n' ...
             '   run crosscheck_response() to document the safeguard.\n']);
end

%% ------------------------------------------------------------------------
function summary = summarize_results(results)
% Format the aggregated statistics as plain text (no raw data leave the
% pipeline beyond the residual pattern indicators).
    lines = {};
    lines{end+1} = sprintf('Number of candidate models fitted: %d', numel(results.models));
    for k = 1:numel(results.models)
        lines{end+1} = sprintf('  Model: %s | AICc = %.2f | R^2 = %.4f', ...
            results.models{k}, results.AICc(k), results.R2(k));
    end
    [best_aicc, ibest] = min(results.AICc);
    lines{end+1} = sprintf('Best model by AICc: %s', results.models{ibest});
    for k = 1:numel(results.models)
        lines{end+1} = sprintf('  delta_AICc(%s) = %.2f', results.models{k}, results.AICc(k) - best_aicc);
    end
    lines{end+1} = 'Parameter estimates (95% CI) of the best model:';
    for k = 1:numel(results.param_names)
        lines{end+1} = sprintf('  %s = %.4g  [%.4g, %.4g]', ...
            results.param_names{k}, results.params(k), results.ci(k,1), results.ci(k,2));
    end
    % Residual pattern indicators (sign of the mean residual in the
    % lowest/middle/highest substrate thirds)
    n = numel(results.residuals);
    t = max(1, floor(n/3));
    lines{end+1} = sprintf('Residual sign pattern (low/mid/high [S]): %+d / %+d / %+d', ...
        sign(mean(results.residuals(1:t))), ...
        sign(mean(results.residuals(t+1:2*t))), ...
        sign(mean(results.residuals(2*t+1:end))));
    if isfield(results, 'v') && isfield(results, 'S')
        lines{end+1} = sprintf('Velocity trend at high [S]: %s', ...
            ternary(results.v(end) < results.v(round(end/2)), 'decreasing', 'non-decreasing'));
    end
    summary = strjoin(lines, newline);
end

%% ------------------------------------------------------------------------
function flags = heuristic_rules(results)
% Hard-coded safeguard rules. Each flag is a falsifiable hypothesis.
    flags = {};
    [~, ibest] = min(results.AICc);
    dAICc = results.AICc - results.AICc(ibest);

    % Rule 1: competing models statistically indistinguishable
    if any(dAICc > 0 & dAICc < 2)
        flags{end+1} = ['Two or more models differ by delta_AICc < 2: they are ' ...
            'statistically indistinguishable; prefer the simpler one or collect ' ...
            'more data before assigning a mechanism.'];
    end
    % Rule 2: wide confidence intervals (parameter identifiability)
    rel_ci = (results.ci(:,2) - results.ci(:,1)) ./ max(abs(results.params), eps);
    if any(rel_ci > 1.0)
        flags{end+1} = ['At least one parameter has a 95% CI wider than its own ' ...
            'estimate: possible parameter identifiability problem; consider ' ...
            'additional substrate levels or a simpler model.'];
    end
    % Rule 3: good R^2 but structured residuals (misfit signature)
    n = numel(results.residuals); t = max(1, floor(n/3));
    signs = [sign(mean(results.residuals(1:t))), ...
             sign(mean(results.residuals(t+1:2*t))), ...
             sign(mean(results.residuals(2*t+1:end)))];
    if max(results.R2) > 0.95 && numel(unique(signs)) > 1
        flags{end+1} = ['R^2 is high but the residual signs change across the ' ...
            'substrate range: systematic misfit despite a good global fit. ' ...
            'Inspect the curve shape before accepting the model.'];
    end
    % Rule 4: velocity decreases at high [S] (substrate inhibition signature)
    if isfield(results, 'v') && isfield(results, 'S')
        if results.v(end) < results.v(round(end/2))
            flags{end+1} = ['The velocity decreases at high substrate ' ...
                'concentrations: the signature of substrate inhibition. ' ...
                'Consider refitting with the Haldane model and re-ranking by AICc.'];
        end
    end
end

%% ------------------------------------------------------------------------
function prompt = build_prompt(summary, flags, stage)
    stage_tasks = { ...
        'Inspect these kinetic results and flag anything unusual in the data quality or experimental design.', ...
        'Based on the curve-shape indicators, propose which kinetic models are plausible candidates and why.', ...
        'Diagnose the fit: interpret the AICc differences, the confidence intervals, and the residual pattern. Recommend whether the current best model is acceptable or which alternative should be refitted.', ...
        'Given this kinetic characterization, suggest reasonable parameter ranges and design considerations for a reactor-level simulation.'};

    lines = {};
    lines{end+1} = 'You are assisting with an enzyme kinetics analysis. Below is a structured summary of numerical results already computed by a deterministic pipeline. Do not invent numerical values; work only with what is provided. Present every recommendation as a hypothesis that the user must verify by rerunning the analysis.';
    lines{end+1} = '';
    lines{end+1} = 'NUMERICAL SUMMARY:';
    lines{end+1} = summary;
    lines{end+1} = '';
    if ~isempty(flags)
        lines{end+1} = 'RULE-BASED WARNINGS ALREADY TRIGGERED (do not simply repeat them; confirm, refine, or refute them with arguments):';
        for k = 1:numel(flags)
            lines{end+1} = sprintf('  - %s', flags{k});
        end
        lines{end+1} = '';
    end
    lines{end+1} = sprintf('YOUR TASK: %s', stage_tasks{stage});
    prompt = strjoin(lines, newline);
end

%% ------------------------------------------------------------------------
function crosscheck_response(llm_keywords, flags)
% Document the safeguard: report which heuristic flags are supported by
% keywords found in the LLM response (e.g., "Haldane", "substrate
% inhibition", "identifiability"). The final decision always belongs to
% the AICc ranking and the user.
    fprintf('--- SAFEGUARD CROSS-CHECK ---\n');
    if isempty(flags)
        fprintf('   No heuristic flags to check against.\n');
        return;
    end
    for k = 1:numel(flags)
        hit = any(cellfun(@(kw) contains(lower(flags{k}), lower(kw)), llm_keywords));
        if hit
            fprintf('   Flag %d: SUPPORTED by the LLM response.\n', k);
        else
            fprintf('   Flag %d: not addressed by the LLM response - verify manually.\n', k);
        end
    end
    fprintf(['Remember: LLM suggestions are hypotheses. The AICc ranking\n' ...
             'and the user make the final decision.\n']);
end

%% ------------------------------------------------------------------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
