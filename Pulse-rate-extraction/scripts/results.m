%% === Remote PPG Analysis: Final Script for Submission ===

clear; close all; clc;

%% === Path Setup ===
project_root = fileparts(mfilename('fullpath'));
cd(project_root);

addpath(genpath(fullfile(project_root, 'routine')));
addpath(genpath(fullfile(project_root, 'data')));

results_dir = fullfile(project_root, 'results');
subdirs = dir(results_dir);
subdirs = subdirs([subdirs.isdir] & ~ismember({subdirs.name}, {'.', '..'}));

% Keep only numeric folders (patient IDs)
is_numeric = cellfun(@(x) ~isempty(regexp(x, '^\d+$', 'once')), {subdirs.name});
subdirs = subdirs(is_numeric);
[~, idx] = sort(str2double({subdirs.name}));
subdirs = subdirs(idx);

% Load patient structs
total_patients = numel(subdirs);
patients = cell(1, total_patients);
for i = 1:total_patients
    folder_path = fullfile(results_dir, subdirs(i).name);
    mat_files = dir(fullfile(folder_path, '*.mat'));

    if isempty(mat_files)
        warning('No .mat file found in folder %s', subdirs(i).name);
        continue;
    end

    s = load(fullfile(folder_path, mat_files(1).name));
    patients{i} = s;
end

%% === Compute Mean SNR Across Methods ===
method_names = {'bvp', 'S_RoverG', 'S_XoverY_b', 'S_XoverY_f', 'S_XsminaYs', 'S_ICA', 'S_PCA'};
num_methods = length(method_names);
SNR_total = zeros(num_methods, 1);

for i = 1:total_patients
    SNR_struct = patients{i}.results.SNR;
    for j = 1:num_methods
        method = method_names{j};
        SNR_value = SNR_struct.(method);
        SNR_total(j) = SNR_total(j) + SNR_value;
    end
end

SNR_mean = SNR_total / total_patients;

% Plot and Save Mean SNR
figure;
bar(SNR_mean(2:end));
set(gca, 'XTickLabel', method_names(2:end), 'XTick', 1:num_methods-1);
ylabel('Mean SNR [dB]');
title('Mean SNR per Method');
grid on;

output_folder = fullfile(results_dir, 'overall_results');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
saveas(gcf, fullfile(output_folder, 'mean_SNR_per_method.png'));

%% === Compute Agreement Metrics ===
methods = {'bvp','RoverG','XoverY_basic','XoverY_fixed','XsminaYs','ICA','PCA'};
M = numel(methods);
pulse_rate = [];
pr_ref_all = [];

start_idx = 1;
for i = 1:total_patients
    N = length(patients{i}.results.HR_reference);
    pr_ref_all = [pr_ref_all, patients{i}.results.HR_reference(:)'];
    for m = 1:M
        actual_length = size(patients{i}.results.HR_estimated, 2);
        if actual_length == N
            pulse_rate(m, start_idx:start_idx+N-1) = patients{i}.results.HR_estimated(m, :);
        end
    end
    start_idx = start_idx + N;
end

% Correlation, regression, std dev and RMSE
r = zeros(1, M);
B = zeros(1, M);
sigma = zeros(1, M);
epsilon = zeros(1, M);

for i = 1:M
    x = pulse_rate(i,:);
    y = pr_ref_all;
    r(i) = corr(x', y');
    coeffs = polyfit(x, y, 1);
    B(i) = coeffs(1);
    diff = x - y;
    sigma(i) = std(diff);
    epsilon(i) = sqrt(mean(diff.^2));
end

%% === Bland-Altman Plots ===
bias = zeros(1, M);
loaU = zeros(1, M);
loaL = zeros(1, M);
outMat = false(M, length(pr_ref_all));
perc_exact = zeros(1, M);

for i = 1:M
    diff = pulse_rate(i,:) - pr_ref_all;
    meanHR = (pulse_rate(i,:) + pr_ref_all)/2;

    bias(i) = mean(diff);
    sd = std(diff);
    loaU(i) = bias(i) + 1.96 * sd;
    loaL(i) = bias(i) - 1.96 * sd;

    outMat(i,:) = (diff < loaL(i)) | (diff > loaU(i));

    % Bland-Altman plot
    figure;
    scatter(meanHR, diff, 'b');
    hold on;
    yline(bias(i), 'r', 'Bias');
    yline([loaL(i), loaU(i)], 'g--', '+/-1.96 SD');
    xlabel('Mean HR [BPM]');
    ylabel('Difference [BPM]');
    title(sprintf('Bland-Altman – %s', methods{i}));
    grid on;
    saveas(gcf, fullfile(output_folder, sprintf('BlandAltman_%s.png', methods{i})));

    delta = pr_ref_all - pulse_rate(i,:);
    perc_exact(i) = sum(abs(delta) < 5) / length(delta) * 100;
    within_limits.(methods{i}) = sum(~outMat(i,:)) / length(delta) * 100;
end

% disp(within_limits);
%% === Outlier Analysis per Subject and Method ===
HRref = [];
HRalg = [];
Subj = [];
Meth = [];

for p = 1:total_patients
    pr_ref = patients{p}.results.HR_reference(:)';
    pr_alg = patients{p}.results.HR_estimated;
    N = numel(pr_ref);
    HRref = [HRref, pr_ref];
    HRalg = [HRalg, pr_alg];
    Subj = [Subj, p * ones(1, N)];
    Meth = [Meth, repmat((1:M)', 1, N)];
end

out_by_subj_method = zeros(total_patients, M);
for subj = 1:total_patients
    for m = 1:M
        idx = (Subj == subj);
        out_flags = outMat(m, idx);
        out_by_subj_method(subj, m) = 100 * mean(out_flags);
    end
end

T = array2table(out_by_subj_method, 'VariableNames', methods, ...
    'RowNames', arrayfun(@(x) sprintf('Subject_%02d', x+10), 1:total_patients, 'UniformOutput', false));

% Save the outlier heatmap
figure('Color','w');
imagesc(out_by_subj_method);
colormap(jet);
colorbar;
caxis([0 max(out_by_subj_method(:))]);
xticks(1:M);
xticklabels(methods);
yticks(1:total_patients);
yticklabels(arrayfun(@(x) sprintf('Subject %02d', x+10), 1:total_patients, 'UniformOutput', false));
xlabel('Method');
ylabel('Subject');
title('Outlier Rate per Subject and Method');
grid on;
box on;

saveas(gcf, fullfile(output_folder, 'outlier_rate_heatmap.png'));

% Save summary statistics
% save(fullfile(output_folder, 'summary_metrics.mat'), 'r', 'B', 'sigma', 'epsilon', 'bias', 'loaU', 'loaL', 'perc_exact', 'within_limits');
