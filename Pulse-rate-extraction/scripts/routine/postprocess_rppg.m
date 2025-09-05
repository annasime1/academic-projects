function [bvp_f, S_RoverG_f, S_XoverY_basic_f, S_XoverY_fixed_f, ...
    S_XsminaYs_f, S_ICA_f, S_PCA_f] = postprocess_rppg( ...
    Fs_target, bvp_sync, S_RoverG, S_XoverY_basic, ...
    S_XoverY_fixed, S_XsminaYs, S_ICA, S_PCA, results_dir, show_plot)

% POSTPROCESS_RPPG: Filters the extracted rPPG signals and generates comparison plots.
%
% INPUT:
%   Fs_target         - Target sampling frequency
%   bvp_sync          - Synchronized BVP signal
%   S_RoverG          - rPPG signal (R over G method)
%   S_XoverY_basic    - rPPG signal (basic X over Y method)
%   S_XoverY_fixed    - rPPG signal (standardized X over Y method)
%   S_XsminaYs        - rPPG signal (Xs - alpha*Ys method)
%   S_ICA             - ICA-derived signal
%   S_PCA             - PCA-derived signal
%   results_dir       - (optional) directory where the plot will be saved
%   show_plot         - (optional) true to display comparison plot
%
% OUTPUT:
%   *_f - Filtered versions of the input signals

if nargin < 10
    show_plot = false;
end

%% === Bandpass filter parameters (40–240 BPM) ===
low_cutoff = 40 / 60;
high_cutoff = 240 / 60;
F_n = Fs_target / 2;

[b, a] = butter(4, [low_cutoff, high_cutoff] / F_n, 'bandpass');

%% === Apply bandpass filter ===
bvp_f            = filtfilt(b, a, bvp_sync);
S_RoverG_f       = filtfilt(b, a, S_RoverG);
S_XoverY_basic_f = filtfilt(b, a, S_XoverY_basic);
S_XoverY_fixed_f = filtfilt(b, a, S_XoverY_fixed);
S_XsminaYs_f     = filtfilt(b, a, S_XsminaYs);
S_ICA_f          = filtfilt(b, a, S_ICA');
S_PCA_f          = filtfilt(b, a, S_PCA);

%% === Ensure row vectors for consistency ===
signals_raw  = {S_RoverG, S_XoverY_basic, S_XoverY_fixed, S_XsminaYs, S_ICA, S_PCA};
signals_filt = {S_RoverG_f, S_XoverY_basic_f, S_XoverY_fixed_f, S_XsminaYs_f, S_ICA_f, S_PCA_f};
titles       = {'R over G', 'X over Y (basic)', 'X over Y (fixed)', 'Xs - \alphaYs', 'ICA', 'PCA'};

samples = 1:min(1150, length(S_RoverG));

%% === Plot original vs filtered signals ===
if show_plot
    figure;
    for i = 1:6
        subplot(3, 2, i);
        plot(samples, signals_raw{i}(samples), 'k'); hold on;
        plot(samples, signals_filt{i}(samples), 'LineWidth', 1.2);
        title(['Signal ', titles{i}, ' - Original vs Filtered']);
        legend('Original', 'Filtered');
        grid on;
    end
    sgtitle('Comparison of Original vs Filtered rPPG Signals');
    %% === Save figure if requested ===
    if nargin >= 9 && isfolder(results_dir)
        saveas(gcf, fullfile(results_dir, 'comparison_filtered_signals.png'));
    end
end


end
