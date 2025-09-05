function [bvp_aligned, signals_aligned, optimal_shift] = align_signals(bvp_sync, rppg_signal, F_s, show_plot, results_path)
% ALIGN_SIGNALS Aligns the rPPG signal with the reference BVP and optionally shows/saves the comparison plot
%
% Input:
%   bvp_sync      - synchronized and filtered BVP signal
%   rppg_signal   - filtered rPPG signal to be aligned
%   F_s           - sampling frequency (Hz)
%   show_plot     - (optional, bool) true to show comparison plot
%   results_path  - (optional, string) folder path to save the figure
%
% Output:
%   bvp_aligned        - BVP trimmed to the same length
%   signals_aligned    - aligned and trimmed rPPG
%   optimal_shift      - estimated delay in number of samples

    if nargin < 4
        show_plot = false;
    end

    if nargin < 5
        results_path = '';
    end

    %% 1. Remove edges of rPPG (caused by filters)
    buffer = 500;
    rppg_signal = rppg_signal(buffer:end - buffer);

    %% 2. Normalization with moving window (~5 seconds)
    w = 5 * F_s;
    overlap = round(w / 2);
    L = length(rppg_signal);
    l = length(bvp_sync);

    S_norm = zeros(1, L);
    bvp_norm = zeros(1, l);

    for i = overlap + 1:L - overlap - 1
        S_norm(i) = (rppg_signal(i) - mean(rppg_signal(i - overlap:i + overlap))) / std(rppg_signal(i - overlap:i + overlap));
        if i <= l - overlap -1
            bvp_norm(i) = (bvp_sync(i) - mean(bvp_sync(i - overlap:i + overlap))) / std(bvp_sync(i - overlap:i + overlap));
        end
    end

    % Normalize edges
    S_norm(1:overlap) = (rppg_signal(1:overlap) - mean(rppg_signal(1:w))) / std(rppg_signal(1:w));
    S_norm(end - overlap:end) = (rppg_signal(end - overlap:end) - mean(rppg_signal(end - w:end))) / std(rppg_signal(end - w:end));
    bvp_norm(1:overlap) = (bvp_sync(1:overlap) - mean(bvp_sync(1:w))) / std(bvp_sync(1:w));
    bvp_norm(end - overlap:end) = (bvp_sync(end - overlap:end) - mean(bvp_sync(end - w:end))) / std(bvp_sync(end - w:end));

    %% 3. Cross-correlation
    [c, lags] = xcorr(S_norm, bvp_norm);
    [~, max_idx] = max(c);
    optimal_shift = abs(lags(max_idx));

    %% 4. Alignment
    if optimal_shift >= 0 
        if optimal_shift + length(bvp_norm) <= length(S_norm)
            signals_aligned = S_norm(optimal_shift + 1:optimal_shift + length(bvp_norm));
            bvp_aligned = bvp_norm(1:length(signals_aligned));
        else
            max_len = length(S_norm) - optimal_shift;
            signals_aligned = S_norm(optimal_shift + 1 : end);
            bvp_aligned = bvp_norm(1 : max_len);
            if (length(S_norm) < L * 0.95)
                disp('Alignment error 1');
            end
        end
    else 
        disp('Alignment error 2');
        optimal_shift = 1;
        signals_aligned = S_norm(1:length(bvp_norm));
        bvp_aligned = bvp_norm(1:length(signals_aligned));
    end

    %% 5. Plot comparison of 1 minute segment
    if show_plot
        start_idx = 30000;                % Start point (sample)
        window_duration = 50;             % Segment duration in seconds
        window_samples = F_s * window_duration;
        end_idx = min(start_idx + window_samples - 1, length(bvp_aligned));

        if end_idx > start_idx  % Ensure valid interval
            t = (start_idx:end_idx) / F_s;

            fig = figure('Name', 'BVP vs rPPG (Aligned)', 'NumberTitle', 'off');
            plot(t, signals_aligned(start_idx:end_idx), 'b', 'DisplayName', 'rPPG (aligned)'); hold on;
            plot(t, bvp_aligned(start_idx:end_idx), 'r', 'DisplayName', 'BVP');
            title(sprintf('Comparison BVP vs rPPG (aligned)'));
            xlabel('Time [s]');
            xlim([min(t), max(t)]);
            ylabel('Amplitude (normalized)');
            legend(); grid on;

            % Save figure
            if isfolder(results_path)
                saveas(fig, fullfile(results_path, 'alignment_bvp_vs_rppg.png'));
            end
        else
            warning('Segment too short for plotting: start = %d, end = %d', start_idx, end_idx);
        end
    end

end
