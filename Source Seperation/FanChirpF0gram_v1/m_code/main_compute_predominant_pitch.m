%% main_compute_predominant_pitch
%
% Program to process a database of audio files computing the FChT and the F0gram
% for each file along with an ordered list of pitch candidates for each frame.
% This estimates can be evaluated against f0 labels to asses its
% performance using other scripts provided: 
%                                     - main_visualize_individual_results.m
%                                     - main_visualize_global_results.m

close all
clear all

%% ============= PARAMETERS AND DESIGN =============

if(isunix); dir_slash = '/'; else dir_slash = '\'; end    
% add path to the c_code
addpath(['.' dir_slash '..' dir_slash 'c_code']);
% add path to set_parameters script
addpath(['.' dir_slash '..' dir_slash]);

% set FChT parameters
set_parameters;

% design for efficient fcht computation
[warps cqt f0s accums] = design_fcht(labels_params.fs_database, nfft, fmax, cqt_params, warp_params, f0_params, glogs_params);


%% ============= DATABASE FILES =============

% struct with the name and complete path to each audio file in the database 
db_files = files_in_path(labels_params.database_path, 'wav');
% number of files in the database
num_db_files = length(db_files);


%% ====== PROCESSING ======

tic
for i_file = 1 : num_db_files

    % audio filename without path
    labels_params.audio_file = db_files(i_file).clean_file_name;

    % mat files for the result of the analysis
    f0gram_filename = [labels_params.path_mat labels_params.audio_file '_f0gram.mat'];
    pitch_filename = [labels_params.path_mat labels_params.audio_file '_pitch.mat'];

    % if pitch detection has NOT been done previously (file does not exist)
    if ~(exist(pitch_filename,'file') == 2)

        % audio file reading
        disp(['Processing file ' labels_params.audio_file ' ...'])
        [y, fs] = audioread(db_files(i_file).name);
        if (labels_params.fs_database ~= fs); disp('WARNING: audio fs different from design fs'); end;
        % only one channel (left)
        if (size(y,2) > 1); y = y(:,1); disp('WARNING: only left channel considered'); end
        
        % fcht computation
        [f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec] = compute_fcht(y, fs, nfft, fmax, hop, cqt, warps, f0s, accums, f0_params, glogs_params, labels_params);
        
        % save results
        save(pitch_filename, 'f0s','f0_hyps_indxs','val_f0_hyps','t');
        if f0_params.save_f0gram; save(f0gram_filename, 'f0gram','t','f0s'); end
                
    else 
        disp(['Skipping file ' labels_params.audio_file ' (remove .mat to force computation)'])
    end
    
end
toc 


