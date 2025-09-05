function [R_filtered, G_filtered, B_filtered, F_s, t_resampled] = preprocess_rgb(RawData, timestamp, Fs_target, show_plot, config)
% PREPROCESS_RGB: Corrects jumps, interpolates, and filters RGB channels.
%
% INPUT:
%   RawData     - Nx3 matrix with columns [R, G, B]
%   timestamp   - original time vector
%   Fs_target   - desired resampling frequency (Hz)
%   show_plot   - (bool, optional) true to display comparison plots
%   config      - configuration structure, used to save plots if enabled
%
% OUTPUT:
%   R_filtered, G_filtered, B_filtered - preprocessed RGB channels
%   F_s             - effective sampling frequency (equal to Fs_target)
%   t_resampled     - uniformly resampled time vector

    if nargin < 4
        show_plot = false;
    end

    %% === Extract RGB channels ===
    R_raw = RawData(:,1);
    G_raw = RawData(:,2);
    B_raw = RawData(:,3);

    %% === Jump correction ===
    R = correct_jumps(R_raw);
    G = correct_jumps(G_raw);
    B = correct_jumps(B_raw);

    %% === Interpolation ===
    new_time = timestamp(1):1/Fs_target:timestamp(end);
    timestamp = timestamp(1:length(R));

    R_interp = interp1(timestamp, R, new_time, 'spline');
    G_interp = interp1(timestamp, G, new_time, 'spline');
    B_interp = interp1(timestamp, B, new_time, 'spline');

    t_resampled = new_time;
    F_s = Fs_target;

    %% === Low-pass filtering ===
    f_cut = 10;
    F_n = F_s / 2;
    [b_lp, a_lp] = butter(4, f_cut / F_n);

    R_filtered = filtfilt(b_lp, a_lp, R_interp);
    G_filtered = filtfilt(b_lp, a_lp, G_interp);
    B_filtered = filtfilt(b_lp, a_lp, B_interp);

    %% === Comparison plot (raw vs filtered) ===
    if show_plot
        N = min(400, length(R_raw));
        figure('Name', 'RGB Channels: Raw vs Preprocessed', 'NumberTitle', 'off');

        subplot(3,1,1);
        plot(R_raw(1:N), 'k'); hold on;
        plot(R_filtered(1:N), 'r');
        title('Red Channel'); ylabel('Intensity'); grid on;
        legend('Original', 'Processed');

        subplot(3,1,2);
        plot(G_raw(1:N), 'k'); hold on;
        plot(G_filtered(1:N), 'g');
        title('Green Channel'); ylabel('Intensity'); grid on;
        legend('Original', 'Processed');

        subplot(3,1,3);
        plot(B_raw(1:N), 'k'); hold on;
        plot(B_filtered(1:N), 'b');
        title('Blue Channel'); ylabel('Intensity'); xlabel('Frame'); grid on;
        legend('Original', 'Processed');

        % === Save figure if results folder exists ===
        if isfolder(config.results_path)
            fig = gcf;  % Get current figure
            saveas(fig, fullfile(config.results_path, 'preprocessing_comparison.png'));
        end
    end
end

function signal_out = correct_jumps(signal_in)
% CORRECT_JUMPS: Removes large discontinuities in a signal
% by detecting and compensating sudden jumps.

    threshold = 3 * std(diff(signal_in));
    jumps = find(abs(diff(signal_in)) > threshold);
    signal_out = signal_in;
    cumulative_shift = 0;

    for i = 1:length(jumps)
        idx = jumps(i) + 1;
        jump_size = signal_out(idx) - signal_out(idx - 1);
        cumulative_shift = cumulative_shift + jump_size;
        signal_out(idx:end) = signal_out(idx:end) - jump_size;
    end
end
