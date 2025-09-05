function [bvp_sync, green_corrected, RawData, timestamp] = load_patient_data(config)
% LOAD_PATIENT_DATA: Loads BVP, corrected green channel, RGB data, and timestamps for a given patient.
%
% INPUT:
%   config - struct containing configuration parameters, including:
%            - patient_id
%            - data_path
%
% OUTPUT:
%   bvp_sync        - synchronized BVP signal (mean-centered)
%   green_corrected - corrected green channel signal
%   RawData         - raw RGB matrix (Nx3)
%   timestamp       - timestamp vector for the video frames

    pid = config.patient_id;

    % === Load BVP signal ===
    bvp_file = fullfile(config.data_path, 'bvp_p1.csv');
    bvp = load(bvp_file);
    bvp_sync = bvp(2:end, pid - 10) - mean(bvp(2:end, pid - 10));

    % === Load corrected green channel ===
    green_file = fullfile(config.data_path, 'green_corrected_p1.csv');
    green_corrected = load(green_file);

    % === Load RGB signal ===
    rgb_folder = fullfile(config.data_path, 'RGB');
    rgb_pattern = sprintf('rgb_raw_%d_*.csv', pid);
    rgb_files = dir(fullfile(rgb_folder, rgb_pattern));
    if isempty(rgb_files)
        error('RGB file not found for patient %d', pid);
    end
    RawData = readmatrix(fullfile(rgb_folder, rgb_files(1).name));

    % === Load timestamp ===
    ts_folder = fullfile(config.data_path, 'Timestamp_p1');
    ts_pattern = sprintf('video_%d_*.csv', pid);
    ts_files = dir(fullfile(ts_folder, ts_pattern));
    if isempty(ts_files)
        error('Timestamp file not found for patient %d', pid);
    end
    timestamp = readmatrix(fullfile(ts_folder, ts_files(1).name), 'NumHeaderLines', 1);
end
