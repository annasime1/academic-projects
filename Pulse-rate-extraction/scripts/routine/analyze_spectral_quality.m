function [SNR_results, pulse_rate] = analyze_spectral_quality(signals_aligned, F_s, show_plot, results_path)
% ANALYZE_SPECTRAL_QUALITY: Computes SNR and estimated HR via spectral peak analysis.
%
% Inputs:
%   signals_aligned - struct with aligned rPPG signals (must include field 'bvp' as reference)
%   F_s             - sampling frequency (Hz)
%   show_plot       - (optional) boolean, true to display and save plot (default: false)
%   results_path    - (optional) folder path to save plots
%
% Outputs:
%   SNR_results     - struct with average SNR values per method
%   pulse_rate      - matrix [num_methods x num_samples] with estimated heart rates (BPM)

    if nargin < 3
        show_plot = false;
    end
    if nargin < 4
        results_path = '';
    end

    methods = fieldnames(signals_aligned);
    num_methods = length(methods);

    bvp = signals_aligned.bvp;
    L = length(bvp);

    x = 256;
    w = round(x / 20 * F_s);          % window length
    overlap = round(w / 2);
    step = 1;
    window = hann(w)';
    Nfft = w;

    % Frequency axis in BPM
    frequencies = (0:Nfft/2-1) * (F_s / Nfft) * 60;
    f_min = 40;
    f_max = 240;
    valid_idx = (frequencies >= f_min) & (frequencies <= f_max);
    start_idx = find(abs(frequencies - f_min) == min(abs(frequencies - f_min)), 1);
    end_idx   = find(abs(frequencies - f_max) == min(abs(frequencies - f_max)), 1);

    % Peak widths for template (number of bins)
    w1 = round(5 / 512 * x / 20 * F_s);   % fundamental peak width
    w2 = round(10 / 512 * x / 20 * F_s);  % harmonic peak width

    % Initialize SNR accumulators and pulse_rate matrix
    SNR = zeros(1, num_methods);
    pulse_rate = zeros(num_methods, L);

    counter = 0;

    for i = overlap:step:L - overlap
        counter = counter + 1;

        % === Reference BVP peak detection ===
        sgn = bvp(i - overlap + 1:i + overlap);
        sgn = sgn(:)' .* window;
        spectrum = abs(fft(sgn, Nfft));
        spectrum = spectrum(1:Nfft/2);

        [~, peak_idx] = max(spectrum(valid_idx));
        fundamental_freq = frequencies(valid_idx);
        fundamental_freq = fundamental_freq(peak_idx);
        harmonic_freq = 2 * fundamental_freq;

        % Find indexes of fundamental and harmonic frequencies
        fundamental_idx = find(abs(frequencies - fundamental_freq) == min(abs(frequencies - fundamental_freq)), 1);
        harmonic_idx = find(abs(frequencies - harmonic_freq) == min(abs(frequencies - harmonic_freq)), 1);

        % Create spectral template (1 at fundamental and harmonic peaks)
        template = zeros(size(frequencies));
        for j = 1:length(frequencies)
            if abs(j - fundamental_idx) <= w1 / 2
                template(j) = 1;
            elseif abs(j - harmonic_idx) <= w2 / 2
                template(j) = 1;
            end
        end

        % Compute SNR for BVP
        NUM = sum((spectrum(start_idx:end_idx) .* template(start_idx:end_idx)).^2);
        DEN = sum((spectrum(start_idx:end_idx) .* (1 - template(start_idx:end_idx))).^2);
        SNR(1) = SNR(1) + 10 * log10(NUM / DEN);
        pulse_rate(1, i) = fundamental_freq;

        % === Compute SNR and HR for other rPPG methods ===
        for m = 2:num_methods
            signal = signals_aligned.(methods{m});
            segment = signal(i - overlap + 1:i + overlap);
            segment = segment(:)' .* window;

            spectrum = abs(fft(segment, Nfft));
            spectrum = spectrum(1:Nfft/2);

            NUM = sum((spectrum(start_idx:end_idx) .* template(start_idx:end_idx)).^2);
            DEN = sum((spectrum(start_idx:end_idx) .* (1 - template(start_idx:end_idx))).^2);
            SNR(m) = SNR(m) + 10 * log10(NUM / DEN);

            % Estimate HR
            [~, peak_idx] = max(spectrum(valid_idx));
            fundamental_freq = frequencies(valid_idx);
            fundamental_freq = fundamental_freq(peak_idx);
            pulse_rate(m, i) = fundamental_freq;
        end
    end

    % Average SNR over windows
    SNR = SNR / counter;

    % Prepare output struct with SNR values
    SNR_results = struct();
    for m = 1:num_methods
        SNR_results.(methods{m}) = SNR(m);
    end

    % Optional plotting
    if show_plot
        figure;
        bar(SNR, 'FaceColor', [0.2 0.4 0.8]);
        xticks(1:num_methods);
        xticklabels(methods);
        ylabel('SNR (dB)');
        title('Spectral SNR of rPPG methods');
        grid on;

        if isfolder(results_path)
            saveas(gcf, fullfile(results_path, 'SNR_comparison.png'));
        end
    end
end
