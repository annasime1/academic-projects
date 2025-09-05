clear all; close all; clc;

% === PATH SETUP ===
project_root = fileparts(mfilename('fullpath'));
cd(project_root);
addpath(genpath(fullfile(project_root, 'routine')));
addpath(genpath(fullfile(project_root, 'data')));
results_root = fullfile(project_root, 'results');

% === CONFIGURATION ===
patient_ids = 11:12;  % List of patients to process
Fs_target = 115;      % Target sampling frequency
F_n = Fs_target / 2;  % Nyquist frequency

% Frequency cutoffs in Hz (40 BPM to 240 BPM)
low_cutoff = 40 / 60;
high_cutoff = 240 / 60;

% Butterworth bandpass filter design (4th order)
[b, a] = butter(4, [low_cutoff high_cutoff] / F_n, 'bandpass');

config = struct();
config.Fs_target       = 115;                              % Target sampling frequency
config.plot_enabled    = false;                             % Enable/disable plotting
config.save_results    = false;                            % Enable/disable result saving
config.data_path       = fullfile(project_root, 'data');
config.routine_path    = fullfile(project_root, 'routine');

error_global = zeros(length(patient_ids), 16);

for p = 1:length(patient_ids)

    % === CONFIGURATION ===

    config.patient_id      = patient_ids(p);                                % Patient ID

    fprintf('\n===== Processing patient ID: %d =====\n', config.patient_id);
    
    % Load patient data
    fprintf('Loading data...\n');
    [bvp_sync, ~, RawData, timestamp] = load_patient_data(config);
    
    % Preprocess RGB signals
    fprintf('Preprocessing RGB signals...\n');
    [R_filtered, G_filtered, B_filtered, F_s, ~] = preprocess_rgb(RawData, timestamp, Fs_target, false, struct());
    L = length(R_filtered);
    
    % Normalize signals with sliding window (~1 second window)
    fprintf('Normalizing RGB signals with sliding window...\n');
    w = round(32/20 * Fs_target);  % window length in samples
    R_norm = zeros(1, L);
    G_norm = zeros(1, L);
    B_norm = zeros(1, L);

    % Center part normalization
    for i = w/2+1 : L-w/2
        R_norm(i) = R_filtered(i) / mean(R_filtered(i - w/2 : i + w/2));
        G_norm(i) = G_filtered(i) / mean(G_filtered(i - w/2 : i + w/2));
        B_norm(i) = B_filtered(i) / mean(B_filtered(i - w/2 : i + w/2));
    end
    
    % Edge normalization
    R_norm(1:w/2) = R_filtered(1:w/2) / mean(R_filtered(1:w));
    G_norm(1:w/2) = G_filtered(1:w/2) / mean(G_filtered(1:w));
    B_norm(1:w/2) = B_filtered(1:w/2) / mean(B_filtered(1:w));

    R_norm(L-w/2+1:L) = R_filtered(L-w/2+1:L) / mean(R_filtered(L-w:L));
    G_norm(L-w/2+1:L) = G_filtered(L-w/2+1:L) / mean(G_filtered(L-w:L));
    B_norm(L-w/2+1:L) = B_filtered(L-w/2+1:L) / mean(B_filtered(L-w:L));
    
    % Compute standardized chrominance signals
    fprintf('Computing chrominance signals...\n');
    X_s = 3 * R_norm - 2 * G_norm;
    Y_s = 1.5 * R_norm + G_norm - 1.5 * B_norm;
    
    % Bandpass filter the chrominance signals
    fprintf('Filtering chrominance signals...\n');
    order_fir = 200;
    b_fir = fir1(order_fir, [low_cutoff high_cutoff] / F_n, 'bandpass');
    a_fir = 1;
    
    X_f = filtfilt(b_fir, a_fir, X_s);
    Y_f = filtfilt(b_fir, a_fir, Y_s);
    
    % Sliding window to compute S_XsminaYs signal
    fprintf('Computing S_XsminaYs signal...\n');
    S_XsminaYs = zeros(1, L);
    window = hann(w)';
    overlap = w / 2;
    n = floor(L / overlap) - 1;
    
    for i = 1 : overlap : n * overlap
        X_f_w = X_f(i : i + w - 1) .* window;
        Y_f_w = Y_f(i : i + w - 1) .* window;
        alpha = std(X_f_w) / std(Y_f_w);
        S_XsminaYs(i : i + w - 1) = S_XsminaYs(i : i + w - 1) + (X_f_w - alpha * Y_f_w);
    end
    % Fix edges
    S_XsminaYs(1 : overlap) = S_XsminaYs(overlap);
    S_XsminaYs(i : L) = S_XsminaYs(i);
    
    % Filter BVP and S_XsminaYs signals with Butterworth bandpass
    fprintf('Filtering BVP and S_XsminaYs signals...\n');
    bvp_f = filtfilt(b, a, bvp_sync);
    S_XsminaYs_f = filtfilt(b, a, S_XsminaYs);
    
    % Synchronize BVP and rPPG signals (cross-correlation)
    fprintf('Aligning signals...\n');
    S1 = S_XsminaYs_f(500:end-500);
    bvp_trimmed = bvp_f;
    
    % Normalize on sliding window (~5 sec)
    w_align = 5 * Fs_target;
    overlap_align = round(w_align / 2);
    
    L_align = length(S1);
    bvp_allignment = zeros(1, L_align);
    S_allignment = zeros(1, L_align);
    
    for i = overlap_align + 1 : L_align - overlap_align - 1
        if i <= length(bvp_trimmed) - overlap_align - 1
            bvp_allignment(i) = (bvp_trimmed(i) - mean(bvp_trimmed(i - overlap_align : i + overlap_align))) / std(bvp_trimmed(i - overlap_align : i + overlap_align));
        end
        S_allignment(i) = (S1(i) - mean(S1(i - overlap_align : i + overlap_align))) / std(S1(i - overlap_align : i + overlap_align));
    end
    
    % Edge normalization
    bvp_allignment(1 : overlap_align) = (bvp_trimmed(1 : overlap_align) - mean(bvp_trimmed(1 : w_align))) / std(bvp_trimmed(1 : w_align));
    S_allignment(1 : overlap_align) = (S1(1 : overlap_align) - mean(S1(1 : w_align))) / std(S1(1 : w_align));
    
    bvp_allignment(end - overlap_align : end) = (bvp_trimmed(end - overlap_align : end) - mean(bvp_trimmed(end - w_align : end))) / std(bvp_trimmed(end - w_align : end));
    S_allignment(end - overlap_align : end) = (S1(end - overlap_align : end) - mean(S1(end - w_align : end))) / std(S1(end - w_align : end));
    
    % Cross-correlation to find optimal lag
    [c, lags] = xcorr(S_allignment, bvp_allignment);
    [~, max_idx] = max(c);
    optimal_shift = lags(max_idx);
    
    % Align signals using optimal lag
    if optimal_shift >= 0
        if optimal_shift + length(bvp_allignment) < length(S_allignment)
            S_aligned = S_allignment(optimal_shift + 1 : optimal_shift + length(bvp_allignment));
        else
            S_aligned = S_allignment(optimal_shift + 1 : end);
            bvp_allignment = bvp_allignment(1 : length(S_aligned));
        end
    else
        warning('Negative lag found; using no lag.');
        optimal_shift = 0;
        S_aligned = S_allignment(1 : length(bvp_allignment));
    end
    
    minLen = min(length(bvp_trimmed), length(S_aligned));
    bvp_f = bvp_trimmed(1:minLen);
    S_aligned = S_aligned(1:minLen);  % Also trim S_aligned accordingly
    S_XsminaYs_f = S_XsminaYs_f(500 : end - 500);
    S_XsminaYs_f = S_XsminaYs_f(optimal_shift + 1 : optimal_shift + length(bvp_f));
    
    % Heart rate reference calculation from BVP
    fprintf('Calculating reference heart rate from BVP...\n');
    L_sig = length(bvp_f);
    bvp_norm = zeros(1, L_sig);
    for i = overlap_align + 1 : L_sig - overlap_align
        bvp_norm(i) = (bvp_f(i) - mean(bvp_f(i - overlap_align : i + overlap_align))) / std(bvp_f(i - overlap_align : i + overlap_align));
    end
    
    % Find peaks in normalized BVP to calculate RR intervals
    [pks, locs] = findpeaks(bvp_norm);
    % Filter peaks above mean peak amplitude
    mean_peak = mean(pks);
    valid_idx = pks > mean_peak;
    pks = pks(valid_idx);
    locs = locs(valid_idx);
    
    % Further refine peaks (local maxima)
    refined_pks = [];
    refined_locs = [];
    k = 1;
    for i = 2:length(pks) - 1
        if pks(i) > (pks(i - 1) + pks(i + 1)) / 2.5
            refined_pks(k) = pks(i);
            refined_locs(k) = locs(i);
            k = k + 1;
        end
    end
    
    RR_intervals = diff(refined_locs) / Fs_target; % in seconds
    
    % Interpolate RR intervals to continuous HR reference
    refined_locs = [1 refined_locs length(bvp_norm)];
    RR_intervals = [RR_intervals(1) RR_intervals(1) RR_intervals RR_intervals(end)];
    RR_interp = interp1(refined_locs, RR_intervals, 1:L_sig, 'spline');
    
    % Calculate heart rate reference (BPM) with sliding window
    pr_ref = zeros(1, L_sig);
    for i = overlap_align : L_sig - overlap_align
        pr_ref(i) = 60 / mean(RR_interp(i - overlap_align + 1 : i + overlap_align));
    end
    pr_ref = pr_ref(750 : end - 750);
    
    % === WINDOW OPTIMIZATION ===
    fprintf('Starting window size optimization...\n');
    
    window_sizes = 32 : 32 : 512;
    errors = zeros(1, length(window_sizes));
    
    % Precompute variables for FFT window
    w_fft = 256 / 20 * Fs_target;
    window_fft = hann(w_fft)';
    overlap_fft = w_fft / 2;
    frequencies = (0 : w_fft / 2 - 1) * (Fs_target / w_fft) * 60; % BPM scale
    
    valid_freq_idx = frequencies >= 40 & frequencies <= 240;
    
    for k = 1 : length(window_sizes)
        w_temp = round(window_sizes(k) / 20 * Fs_target);
        overlap_temp = round(w_temp / 2);
        S_temp = zeros(1, L);
        
        max_start = length(X_f) - w_temp + 1;
        for i = 1 : overlap : max_start
            X_win = X_f(i : i + w_temp - 1) .* hann(w_temp)';
            Y_win = Y_f(i : i + w_temp - 1) .* hann(w_temp)';
            alpha = std(X_win) / std(Y_win);
            S_temp(i : i + w_temp - 1) = S_temp(i : i + w_temp - 1) + (X_win - alpha * Y_win);
        end
        S_temp(1 : overlap_temp) = S_temp(overlap_temp);
        S_temp(i : end) = S_temp(i);
        
        S_temp_filt = filtfilt(b, a, S_temp);
        S_temp_filt = S_temp_filt(500 : end - 500);
        S_temp_filt = S_temp_filt(optimal_shift + 1 : optimal_shift + length(bvp_f));
        
        L_temp = length(S_temp_filt);
        pr_temp = zeros(1, L_temp);
        
        % Estimate HR via FFT peak detection for each window in S_temp_filt
        for idx = overlap_fft + 1 : L_temp - overlap_fft
            segment = S_temp_filt(idx - overlap_fft + 1 : idx + overlap_fft) .* window_fft;
            spectrum = abs(fft(segment, w_fft));
            [~, peak_idx] = max(spectrum(valid_freq_idx));
            freqs_in_band = frequencies(valid_freq_idx);  % extract frequencies in the band
            pr_temp(idx) = freqs_in_band(peak_idx);
        end
        
        pr_temp = pr_temp(750 : end - 750);
        errors(k) = sum(abs(pr_ref - pr_temp));
        
        fprintf('Window size %d processed, error = %.4f\n', window_sizes(k), errors(k)/length(pr_ref));
    end
    
    % Find window with minimum error
    [min_error, min_idx] = min(errors);
    best_window = window_sizes(min_idx);
    bw_seconds = best_window / 20;
    
    fprintf('Patient %d - Optimal window size: %d samples (%.2f seconds), min error = %.4f\n', config.patient_id, best_window, bw_seconds, min_error);
    
    % Plot error curve
    figure;
    plot(window_sizes ./20, errors, '-o');
    title(sprintf('Window size optimization error for patient %d', config.patient_id));
    xlabel('Window size (seconds)');
    ylabel('Sum absolute error');
    grid on;
    error_global(p, :) = errors;
end
mae = error_global ./ length(bvp_f);
fprintf('\nProcessing completed for all patients.\n');



