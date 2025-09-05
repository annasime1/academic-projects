function results = run_pipeline(config)
% RUN_PIPELINE Executes the complete rPPG analysis pipeline for a single patient.
%
% INPUT:
%   config - Structure containing the following required fields:
%       config.patient_id      : (int) Patient ID number to analyze.
%       config.Fs_target       : (double) Target sampling frequency (e.g., 115 Hz).
%       config.routine_path    : (string) Path to the routines/functions directory.
%       config.results_path    : (string) Destination path for saving results.
%       config.plot_enabled    : (bool) To enable/disable plot
%
% OUTPUT:
%   results - Structure containing:
%       .patient_id        : Patient ID
%       .signals_rppg      : Struct with filtered and aligned rPPG & BVP signals
%       .HR_reference      : Reference heart rate (from BVP)
%       .HR_estimated      : Estimated HR from all rPPG methods
%       .SNR               : Spectral SNR values per method
%       .agreement         : Statistical agreement metrics vs. reference HR
%

fprintf('=== Starting analysis for patient %d ===\n', config.patient_id);

%% === 1. Load patient data ===
fprintf('[1] Loading patient data...\n');
[bvp_sync, ~ , RawData, timestamp] = load_patient_data(config);

%% === 2. Preprocess RGB channels ===
fprintf('[2] Preprocessing RGB signals...\n');
[R_filt, G_filt, B_filt, F_s, ~ ] = preprocess_rgb(RawData, timestamp, config.Fs_target, config.plot_enabled, config);

%% === 3. Extract rPPG signals ===
fprintf('[3] Extracting rPPG signals using various methods...\n');
[S_RoverG, S_XoverY_basic, S_XoverY_fixed, S_XsminaYs] = compute_rppg_methods(R_filt, G_filt, B_filt, F_s, 32);
[S_ICA, S_PCA] = extract_ica_pca(R_filt, G_filt, B_filt, F_s);

%% === 4. Postprocess rPPG signals ===
fprintf('[4] Postprocessing rPPG and BVP signals...\n');
[bvp_f, S_RoverG_f, S_XoverY_basic_f, S_XoverY_fixed_f, S_XsminaYs_f, S_ICA_f, S_PCA_f] = postprocess_rppg(...
    config.Fs_target, bvp_sync, S_RoverG, S_XoverY_basic, S_XoverY_fixed, S_XsminaYs, S_ICA, S_PCA, config.results_path, config.plot_enabled);

%% === 5. Align rPPG with BVP ===
fprintf('[5] Aligning rPPG signal with BVP reference...\n');
[bvp_aligned , ~ , shift] = align_signals(bvp_f, S_XsminaYs_f, F_s, config.plot_enabled, config.results_path);

%% === 6. Trimming and shifting all rPPG signals ===
fprintf('[6] Trimming and shifting all rPPG signals...\n');
l = length(bvp_aligned);
bvp_f = bvp_f(1:l)';

% Trim 500 samples from the beginning and end
S_RoverG_f       = S_RoverG_f(500:end-500);
S_XoverY_basic_f = S_XoverY_basic_f(500:end-500);
S_XoverY_fixed_f = S_XoverY_fixed_f(500:end-500);
S_XsminaYs_f     = S_XsminaYs_f(500:end-500);
S_ICA_f          = S_ICA_f(500:end-500);
S_PCA_f          = S_PCA_f(500:end-500);

% Align all signals using the computed shift
S_RoverG_f       = S_RoverG_f(shift + 1 : shift + l);
S_XoverY_basic_f = S_XoverY_basic_f(shift + 1 : shift + l);
S_XoverY_fixed_f = S_XoverY_fixed_f(shift + 1 : shift + l);
S_ICA_f          = S_ICA_f(shift + 1 : shift + l);
S_PCA_f          = S_PCA_f(shift + 1 : shift + l);
S_XsminaYs_f     = S_XsminaYs_f(shift + 1 : shift + l);

% Store aligned signals
signals_aligned = struct();
signals_aligned.bvp          = bvp_f;
signals_aligned.S_RoverG     = S_RoverG_f;
signals_aligned.S_XoverY_b   = S_XoverY_basic_f;
signals_aligned.S_XoverY_f   = S_XoverY_fixed_f;
signals_aligned.S_XsminaYs   = S_XsminaYs_f;
signals_aligned.S_ICA        = S_ICA_f;
signals_aligned.S_PCA        = S_PCA_f;

%% === 7. Extract reference heart rate from BVP ===
fprintf('[7] Extracting reference heart rate from BVP...\n');
[pr_ref, ~ ] = extract_reference_hr(bvp_f, F_s, config.plot_enabled);

%% === 8. Compute HR from rPPG ===
fprintf('[8] Computing FFT-based HR estimation...\n');
pulse_rate = compute_hr_fft(signals_aligned, pr_ref, F_s, config.results_path, config.plot_enabled);

%% === 9. Spectral quality analysis ===
fprintf('[9] Performing spectral SNR analysis on rPPG methods...\n');
SNR_results = analyze_spectral_quality(signals_aligned, F_s, config.plot_enabled, config.results_path);

%% === 10. Agreement analysis ===
fprintf('[10] Performing agreement analysis between rPPG and reference HR...\n');
agreement_stats = analyze_agreement(pulse_rate, signals_aligned, config.plot_enabled, config.results_path);

%% === Output structure ===
results = struct();
results.patient_id      = config.patient_id;
results.signals_rppg    = signals_aligned;
results.HR_reference    = pr_ref;
results.HR_estimated    = pulse_rate;
results.SNR             = SNR_results;
results.agreement       = agreement_stats;

fprintf('=== Analysis completed for patient %d ===\n', config.patient_id);
end
