function [S_ICA, S_PCA] = extract_ica_pca(R_filtered, G_filtered, B_filtered, F_s)
% EXTRACT_ICA_PCA Extracts rPPG signals using ICA and PCA on sliding windows
%
% INPUT:
%   R_filtered, G_filtered, B_filtered - filtered RGB channels
%   F_s - sampling frequency (Hz)
%
% OUTPUT:
%   S_ICA - ICA-based extracted rPPG signal
%   S_PCA - PCA-based extracted rPPG signal

%% === ICA - Independent Component Analysis ===

L = length(R_filtered);
w = round(512 / 20 * F_s);      % Window size ~1 second
window = hann(w);               % Hann window
overlap = w / 2;
n = floor(L / overlap) - 1;
num_components = 3;

S_ICA = zeros(L, 1);

% Detrending using Smoothing Priors
T = 1 / F_s;
F_n = F_s / 2;
f_cut = 0.89;
lambda = (1 - cos(2 * pi * f_cut * T))^4 / (16 * (sin(2 * pi * f_cut * T))^2);
D = spdiags(ones(L, 1) * [1 -2 1], [0 1 2], L - 2, L);

R_detr = ((speye(L) + lambda * (D' * D)) \ R_filtered')';
G_detr = ((speye(L) + lambda * (D' * D)) \ G_filtered')';
B_detr = ((speye(L) + lambda * (D' * D)) \ B_filtered')';

RGB_detr = [R_detr' G_detr' B_detr'];

for i = 1:overlap:n * overlap
    RGB_w = RGB_detr(i:i + w - 1, :) .* repmat(window, 1, 3);
    
    mu = mean(RGB_w, 1);
    sigma = std(RGB_w, 0, 1);
    RGB_norm = (RGB_w - mu) ./ sigma;

    [ICs, ~, ~, ~] = fastICA(RGB_norm', num_components, 'kurtosis', 0);
    ICs = ICs';

    % Select component with highest peak in valid heart rate band
    max_peak = 0;
    best_component = 1;
    for j = 1:num_components
        Y = abs(fft(ICs(:, j)));
        Y = Y(1:floor(end / 2));
        freqs = (0:length(Y) - 1) * (F_s / length(ICs(:, j)));
        valid_idx = (freqs >= 0.67) & (freqs <= 4);
        max_Y = max(Y(valid_idx));
        if max_Y > max_peak
            max_peak = max_Y;
            best_component = j;
        end
    end

    ICs(:, best_component) = movmean(ICs(:, best_component), 5);
    len = min(w, L - i + 1);

    if size(ICs, 1) >= len
        S_ICA(i:i + len - 1) = S_ICA(i:i + len - 1) + ICs(1:len, best_component);
    else
        warning('ICA window skipped at index %d (size mismatch)', i);
    end
end

S_ICA = fillmissing(S_ICA, 'linear');
S_ICA(isinf(S_ICA)) = 0;

%% === PCA - Principal Component Analysis ===

S_PCA = zeros(L, 1);
low_cutoff = 40 / 60;
high_cutoff = 240 / 60;

b = fir1(250, [low_cutoff, high_cutoff] / F_n, 'bandpass');
a = 1;

R_f = filtfilt(b, a, R_filtered);
G_f = filtfilt(b, a, G_filtered);
B_f = filtfilt(b, a, B_filtered);

RGB_filt = [R_f' G_f' B_f'];

for i = 1:overlap:L - w
    if i + w - 1 > L, break; end

    RGB_w = RGB_filt(i:i + w - 1, :) .* repmat(window, 1, 3);
    [Z, ~, ~, ~] = PCA(RGB_w', num_components);
    Z = Z';

    % Select component with highest peak in valid heart rate band
    max_peak = 0;
    best_component = 1;
    for j = 1:num_components
        Y = abs(fft(Z(:, j)));
        Y = Y(1:floor(end / 2));
        freqs = (0:length(Y) - 1) * (F_s / length(Z(:, j)));
        valid_idx = (freqs >= 0.67) & (freqs <= 4);
        max_Y = max(Y(valid_idx));
        if max_Y > max_peak
            max_peak = max_Y;
            best_component = j;
        end
    end

    len = min(w, L - i + 1);
    S_PCA(i:i + len - 1) = S_PCA(i:i + len - 1) + Z(1:len, best_component);
end

end
