function [S_RoverG, S_XoverY_basic, S_XoverY_fixed, S_XsminaYs] = compute_rppg_methods(R_filtered, G_filtered, B_filtered, F_s, window_size)
% COMPUTE_RPPG_METHODS Computes rPPG signals using classic and advanced methods
%
% Input:
%   R_filtered, G_filtered, B_filtered - filtered RGB channels (not normalized)
%   F_s - sampling frequency (Hz)
%   window_size - window size in frames (used for normalization)
%
% Output:
%   S_RoverG         - G/R method - basic
%   S_XoverY_basic   - basic X/Y method
%   S_XoverY_fixed   - standardized method with fixed weights
%   S_XsminaYs       - advanced method with dynamic alpha compensation

%% === Sliding Window Normalization (~1s) ===
w = round(window_size / 20 * F_s);  % approx. 1.6 seconds
L = length(R_filtered);

R_norm = zeros(1, L);
G_norm = zeros(1, L);
B_norm = zeros(1, L);

for i = w/2 + 1 : L - w/2
    R_norm(i) = R_filtered(i) / mean(R_filtered(i - w/2 : i + w/2));
    G_norm(i) = G_filtered(i) / mean(G_filtered(i - w/2 : i + w/2));
    B_norm(i) = B_filtered(i) / mean(B_filtered(i - w/2 : i + w/2));
end

R_norm(1:w/2) = R_filtered(1:w/2) / mean(R_filtered(1:w));
G_norm(1:w/2) = G_filtered(1:w/2) / mean(G_filtered(1:w));
B_norm(1:w/2) = B_filtered(1:w/2) / mean(B_filtered(1:w));

R_norm(L - w/2 + 1 : L) = R_filtered(L - w/2 + 1 : L) / mean(R_filtered(L - w : L));
G_norm(L - w/2 + 1 : L) = G_filtered(L - w/2 + 1 : L) / mean(G_filtered(L - w : L));
B_norm(L - w/2 + 1 : L) = B_filtered(L - w/2 + 1 : L) / mean(B_filtered(L - w : L));

%% === 1. RoverG Method ===
S_RoverG = G_norm ./ R_norm - 1;

%% === 2. XoverY Basic Method ===
X = R_filtered - G_filtered;
Y = 0.5 * R_filtered + 0.5 * G_filtered - B_filtered;

X_norm = zeros(1, L);
Y_norm = zeros(1, L);

for i = w/2 + 1 : L - w/2
    X_norm(i) = X(i) / mean(X(i - w/2 : i + w/2));
    Y_norm(i) = Y(i) / mean(Y(i - w/2 : i + w/2));
end

X_norm(1:w/2) = X(1:w/2) / mean(X(1:w));
Y_norm(1:w/2) = Y(1:w/2) / mean(Y(1:w));

X_norm(L - w/2 + 1 : L) = X(L - w/2 + 1 : L) / mean(X(L - w : L));
Y_norm(L - w/2 + 1 : L) = Y(L - w/2 + 1 : L) / mean(Y(L - w : L));

S_XoverY_basic = X_norm ./ Y_norm - 1;
S_XoverY_basic(~isfinite(S_XoverY_basic)) = 0;

%% === 3. XoverY Fixed Method ===
% R_s = 0.7682 * R_norm;
% G_s = 0.5121 * G_norm;
% B_s = 0.3841 * B_norm;

X_s = 3 * R_norm - 2 * G_norm;
Y_s = 1.5 * R_norm + G_norm - 1.5 * B_norm;

S_XoverY_fixed = X_s ./ Y_s - 1;

%% === 4. XsminaYs Method ===
low_cutoff = 40 / 60;
high_cutoff = 240 / 60;
F_n = F_s / 2;

order_fir = 200;
b_fir = fir1(order_fir, [low_cutoff high_cutoff] / F_n, 'bandpass');
a_fir = 1;

X_f = filtfilt(b_fir, a_fir, X_s);
Y_f = filtfilt(b_fir, a_fir, Y_s);

S_XsminaYs = zeros(1, L);
window = hann(w);
overlap = w / 2;
n = floor(L / overlap) - 1;

for i = 1:overlap:n * overlap
    if i + w - 1 > L
        break;
    end

    X_w = X_f(i:i + w - 1) .* window';
    Y_w = Y_f(i:i + w - 1) .* window';

    alpha = std(X_w) / std(Y_w);
    S_XsminaYs(i:i + w - 1) = S_XsminaYs(i:i + w - 1) + (X_w - alpha * Y_w);
end

% Fix signal edges
S_XsminaYs(1:overlap) = S_XsminaYs(overlap);
S_XsminaYs(i:L) = S_XsminaYs(i);

end
