function pulse_rate = compute_hr_fft(signals_struct, pr_ref, F_s, results_path, show_plot)
% COMPUTE_HR_FFT Computes heart rate from rPPG signals using spectral (FFT) analysis
%
% Inputs:
%   signals_struct - struct with rPPG signals from different methods (fields are method names)
%   pr_ref         - reference pulse rate signal (e.g., from BVP or ECG)
%   F_s            - sampling frequency (Hz)
%   results_path   - folder path to save plots (optional)
%   show_plot      - boolean flag to show plots (default: true)
%
% Output:
%   pulse_rate     - matrix of estimated pulse rates [num_methods x signal_length]

if nargin < 5
    show_plot = false;
end

methods = fieldnames(signals_struct);
num_methods = length(methods);
L = length(signals_struct.(methods{1}));

% Window parameters for FFT analysis
w = round(256 / 20 * F_s);      % Window length (~256/20 seconds scaled by Fs)
window = hann(w)';              % Hann window
overlap = round(w / 2);         % 50% overlap
step = 1;                      % Step size for sliding window

% Initialize pulse rate matrix
pulse_rate = zeros(num_methods, L);

% Frequency range for heart rate (in BPM)
freq_min = 40;
freq_max = 240;

for i = overlap:step:L - overlap
    for m = 1:num_methods
        signal = signals_struct.(methods{m});
        segment = signal(i - overlap + 1:i + overlap);
        segment = detrend(segment(:))' .* window;  % Detrend and window
        
        N_fft = w;  
        spectrum = abs(fft(segment, N_fft));
        spectrum = spectrum(1:N_fft/2);
        
        frequencies = (0:N_fft/2-1) * (F_s / N_fft) * 60;  % Convert to BPM
        valid_idx = (frequencies >= freq_min) & (frequencies <= freq_max);
        
        [~, peak_idx] = max(spectrum(valid_idx));
        peak_freq = frequencies(valid_idx);
        fundamental_freq = peak_freq(peak_idx);
        
        pulse_rate(m, i) = fundamental_freq;
    end
end

% Remove edge effects by trimming 750 samples at start and end
if L > 1500
    pulse_rate = pulse_rate(:, 751:end - 750);
end

% === Plot results if requested ===
if show_plot
    % Plot HR from first method vs reference
    figure;
    plot(pr_ref, 'r', 'LineWidth', 1.5); hold on;
    plot(pulse_rate(1, :), 'g', 'LineWidth', 1.5);
    legend('Reference Pulse Rate', methods{1}, 'Location', 'best');
    title('Pulse Rate: Reference vs First rPPG Method');
    xlabel('Samples');
    ylabel('Pulse Rate [BPM]');
    grid on;
    
    if nargin >= 4 && isfolder(results_path)
        saveas(gcf, fullfile(results_path, 'HR_fft_vs_reference.png'));
    end

    % Plot all methods comparison (up to 6 methods)
    figure;
    max_methods_to_plot = min(num_methods, 7);
    method_titles = {'R over G', 'X over Y basic', 'X over Y fixed', ...
        'X - αY', 'ICA', 'PCA'};
    
    for m = 2:max_methods_to_plot
        subplot(3, 2, m - 1);
        plot(pr_ref, 'r', 'LineWidth', 1.2); hold on;
        plot(pulse_rate(m, :), 'g', 'LineWidth', 1.2);
        plot(pulse_rate(1, :), 'k--', 'LineWidth', 1);
        title([method_titles{m - 1} ': HR Comparison']);
        legend('Reference', methods{m}, methods{1}, 'Location', 'best');
        xlabel('Samples'); ylabel('Pulse Rate [BPM]');
        grid on;
    end
    sgtitle('Pulse Rate Comparison Across Different Methods');

    if nargin >= 4 && isfolder(results_path)
        saveas(gcf, fullfile(results_path, 'HR_fft_all_methods.png'));
    end
end

end
