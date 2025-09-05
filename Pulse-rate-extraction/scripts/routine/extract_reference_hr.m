function [pr_ref, peak_locs] = extract_reference_hr(bvp_signal, F_s, show_plot, results_path)
% EXTRACT_REFERENCE_HR Estimates heart rate from the BVP signal
%
% INPUT:
%   bvp_signal   - aligned and filtered BVP signal
%   F_s          - sampling frequency (Hz)
%   show_plot    - (bool) if true, generates and saves plots
%   results_path - (string, optional) path to save the plots
%
% OUTPUT:
%   pr_ref     - interpolated heart rate [BPM]
%   peak_locs  - locations of filtered peaks

    if nargin < 3
        show_plot = false;
    end

    if nargin < 4
        results_path = '';
    end

    %% === 1. Moving window normalization ===
    L = length(bvp_signal);
    w = 5 * F_s;
    overlap = round(w / 2);
    bvp_norm = zeros(1, L);

    for i = overlap + 1:L - overlap
        bvp_norm(i) = (bvp_signal(i) - mean(bvp_signal(i - overlap:i + overlap))) / std(bvp_signal(i - overlap:i + overlap));    
    end

    %% === 2. Peak detection on BVP (minmax) ===
    [peaks, locs] = findpeaks(minmax(bvp_norm));

    %% === 3. Amplitude filtering ===
    pks = [];
    lcs = [];
    k = 1;
    mean_val = mean(minmax(bvp_norm));
    for i = 2:length(peaks) - 1
        if peaks(i) > mean_val
            pks(k) = peaks(i);
            lcs(k) = locs(i);
            k = k + 1;
        end
    end

    %% === 4. Local prominence filtering ===
    peaks = pks;
    locs = lcs;
    pks = [];
    lcs = [];
    k = 1;
    for i = 2:length(peaks) - 1
        if (peaks(i) > (peaks(i - 1) + peaks(i + 1)) / 2.5)
            pks(k) = peaks(i);
            lcs(k) = locs(i);
            k = k + 1;
        end
    end

    peak_locs = lcs;

    %% === 5. Plotting detected peaks (optional) ===
    window_start = 1000;
    window_end = 1000 + 115 * 50;

    valid_peaks_idx = (lcs >= window_start) & (lcs <= window_end);
    lcs_filtered = lcs(valid_peaks_idx);
    pks_filtered = pks(valid_peaks_idx);

    if show_plot
        fig1 = figure('Name', 'BVP Peaks Detection');
        plot(window_start:window_end, minmax(bvp_norm(window_start:window_end)), 'b'); hold on;
        plot(lcs_filtered, pks_filtered, 'ro');
        title('BVP Peak Detection in Selected Window');
        xlabel('Samples');
        ylabel('Amplitude');
        legend('Inverted BVP', 'Refined Peaks');
        grid on;

        if isfolder(results_path)
            saveas(fig1, fullfile(results_path, 'bvp_peaks_detection.png'));
        end
    end

    %% === 6. RR intervals and interpolation ===
    RR = diff(lcs) / F_s;
    lcs_ext = [1 lcs L];
    RR_ext = [RR(1) RR(1) RR RR(end)];
    RR = interp1(lcs_ext, RR_ext, 1:L, 'spline');

    %% === 7. Heart Rate estimation (moving window) ===
    w_hr = round(256 / 20 * F_s);
    overlap_hr = round(w_hr / 2);
    pr_ref = zeros(1, L);

    for i = overlap_hr:L - overlap_hr
        pr_ref(i) = 1 / mean(RR(i - overlap_hr + 1:i + overlap_hr)) * 60;
    end

    % Remove edge samples
    pr_ref = pr_ref(750:end - 750);

    %% === 8. Plotting heart rate ===
    if show_plot
        fig2 = figure('Name', 'Reference Pulse Rate');
        t = (1:length(pr_ref)) / F_s;
        plot(t, pr_ref, 'r', 'LineWidth', 1.5);
        title('Reference Pulse Rate');
        xlabel('Time [s]');
        ylabel('HR [BPM]');
        grid on;
        xlim([t(1), t(end)]);

        if isfolder(results_path)
            saveas(fig2, fullfile(results_path, 'reference_hr_plot.png'));
        end
    end
end
