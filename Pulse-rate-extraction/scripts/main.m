%% ============================
% MAIN.M – Configuration & Execution
% scripts/main.m
%% ============================

clc;
clear;
close all;

% === PATH SETUP ===

% Get the root path of the project
project_root = fileparts(mfilename('fullpath'));
cd(project_root);

% Add subfolders to MATLAB path
addpath(genpath(fullfile(project_root, 'routine')));
addpath(genpath(fullfile(project_root, 'data')));

% === PATIENT LOOP ===
for i = 11:11
    close all; clc;

    % === CONFIGURATION ===
    config = struct();
    config.patient_id      = i;                                % Patient ID
    config.Fs_target       = 115;                              % Target sampling frequency
    config.plot_enabled    = false;                             % Enable/disable plotting
    config.save_results    = false;                            % Enable/disable result saving
    config.data_path       = fullfile(project_root, 'data');
    config.routine_path    = fullfile(project_root, 'routine');

    % === Set up results folder for this patient ===
    results_root = fullfile(project_root, 'results');
    results_path = fullfile(results_root, num2str(config.patient_id));  % e.g., results/11/
    if ~exist(results_path, 'dir')
        mkdir(results_path);
    end
    config.results_path = results_path;

    fprintf("Configuration completed for patient %d\n", config.patient_id);

    %% === RUN PIPELINE ===

    results = run_pipeline(config);

    %% === SAVE RESULTS ===

    if config.save_results
        fname = sprintf('results_patient_%d.mat', config.patient_id);
        save(fullfile(config.results_path, fname), 'results');
    end

    fprintf("End of pipeline for patient %d\n", config.patient_id);
end
