function agreement_stats = analyze_agreement(pulse_rate, signals_aligned, show_plot, results_path)
% ANALYZE_AGREEMENT Compares estimated HR with reference HR (BVP)
% using Bland-Altman analysis and % agreement
%
% Input:
%   pulse_rate      - [num_methods x N] HR matrix (1st row = BVP)
%   signals_aligned - struct with methods (including 'bvp')
%   show_plot       - (bool) whether to show plots
%   results_path    - (string) path to save the plots
%
% Output:
%   agreement_stats - struct with bias, std, LOA, % agreement, etc.

if nargin < 3
    show_plot = false;
end
if nargin < 4
    results_path = '';
end

method_names = fieldnames(signals_aligned);
% ref_name = method_names{1};  % 'bvp'
pr_ref = pulse_rate(1, :);
N = length(pr_ref);

agreement_stats = struct();

for m = 2:length(method_names)
    method = method_names{m};
    pr_est = pulse_rate(m, :);
    delta = pr_ref - pr_est;
    mean_val = (pr_ref + pr_est) / 2;

    bias = mean(delta);
    sd = std(delta);
    loa_upper = bias + 1.96 * sd;
    loa_lower = bias - 1.96 * sd;

    perc_within_LOA = sum(delta >= loa_lower & delta <= loa_upper) / N * 100;
    perc_exact = sum(abs(delta) < 5) / N * 100;

    agreement_stats.(method).bias = bias;
    agreement_stats.(method).std = sd;
    agreement_stats.(method).loa = [loa_lower, loa_upper];
    agreement_stats.(method).perc_within_limits = perc_within_LOA;
    agreement_stats.(method).perc_exact_match = perc_exact;

    % === Bland-Altman plot ===
    if show_plot
        % Identify points inside and outside limits of agreement
        in_loa_idx = delta >= loa_lower & delta <= loa_upper;
        out_loa_idx = ~in_loa_idx;

        figure;

        % Points inside limits (blue, transparent, bigger)
        scatter(mean_val(in_loa_idx), delta(in_loa_idx), 20, 'b', 'filled', ...
            'MarkerFaceAlpha', 0.5);
        hold on;

        % Points outside limits (red, filled, more visible)
        scatter(mean_val(out_loa_idx), delta(out_loa_idx), 30, 'r', 'filled');

        % Lines for bias and limits of agreement
        yline(bias, 'r', 'LineWidth', 2, 'Label', 'Bias');
        yline(loa_upper, 'g--', 'LineWidth', 1.5, 'Label', '+1.96 SD');
        yline(loa_lower, 'g--', 'LineWidth', 1.5, 'Label', '-1.96 SD');

        title(['Bland-Altman - ' method]);
        xlabel('Mean of Measurements');
        ylabel('Difference');
        grid on;

        % Save the figure
        if isfolder(results_path)
            saveas(gcf, fullfile(results_path, ['BA_' method '.png']));
        end
    end

end

% === Output to console ===
disp('=== Agreement Statistics ===');
for m = 2:length(method_names)
    method = method_names{m};
    s = agreement_stats.(method);
    fprintf('%s: Bias = %.2f ± %.2f, In LOA = %.1f%%, Exact < 5 BPM = %.1f%%\n', ...
        method, s.bias, 1.96 * s.std, s.perc_within_limits, s.perc_exact_match);
end
end
