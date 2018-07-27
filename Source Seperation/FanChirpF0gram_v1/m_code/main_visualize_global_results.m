%% main_visualize_global_results
%
% Program to compare pitch candidates with f0 labels for an audio database.
% Performance is computed using a soft thresholding (see paper DAFx 2010).

close all;
clear all;

%% ============= PARAMETERS =============

if(isunix); dir_slash = '/'; else dir_slash = '\'; end    
% add path to the c_code
addpath(['.' dir_slash '..' dir_slash 'c_code']);
% add path to set_parameters script
addpath(['.' dir_slash '..' dir_slash]);

% set FChT parameters
set_parameters;

% filename to save global results
results_filename = [labels_params.path_mat 'global_results.mat'];

% evaluation tolerance parameters
tolerance = 2;
tol_min = 1;
tol_max = 3;
min_valid_f0 = 40;

%% ============= DATABASE FILES =============

% struct with the name and complete path to each audio file in the database 
db_files = files_in_path(labels_params.database_path, 'wav');
% number of files in the database
num_db_files = length(db_files);
% files names as cells to be used in chart axis
names_cell = cell(1,num_db_files);

%% ====== PROCESSING ======
% percentage of hits computation

percentage_of_hits = zeros(num_db_files,f0_params.num_f0_hyps);
total_global = zeros(1,f0_params.num_f0_hyps);
corrects_global = zeros(1,f0_params.num_f0_hyps);

% performance
pitch_melody_global = [];
true_pitch_melody_global = [];
val_melody_global = [];

for i_file = 1 : num_db_files
    
    % audio filename without path
    labels_params.audio_file = db_files(i_file).clean_file_name;
    names_cell(i_file) = {labels_params.audio_file};
    disp(['Processing file ' labels_params.audio_file ' ...'])
    
    % load mat file with pitch candidates
    pitch_filename = [labels_params.path_mat labels_params.audio_file '_pitch.mat'];
    load(pitch_filename);

    % load labels
    [f0_value f0_time] = read_f0_labels(labels_params.database_path, labels_params.audio_file, labels_params.database);

    % interpolation of melody frecuency values at frame time instants
    true_pitch = interp1(f0_time, f0_value, t);
    valid_pitch = interp1(f0_time, f0_value > min_valid_f0, t);

    % percetage of hits
    inds_present = find(valid_pitch == 1);
    total_frames = length(inds_present);
    true_pitch_melody = true_pitch(inds_present);

    pitch_melody = f0s(f0_hyps_indxs(1:f0_params.num_f0_hyps,inds_present));
    true_pitch_melody = repmat(true_pitch_melody,[f0_params.num_f0_hyps,1]);

    pitch_melody_global = [pitch_melody_global pitch_melody];
    true_pitch_melody_global = [true_pitch_melody_global true_pitch_melody];

    val_melody_global = [val_melody_global val_f0_hyps(1:f0_params.num_f0_hyps,inds_present)];

    % soft thresholding evaluation 
    diff_pitch = abs(true_pitch_melody-pitch_melody)./true_pitch_melody*100;
    corrects = min(1, max(0, (diff_pitch - tol_max) ./ (tol_min - tol_max)));
    union_corrects = corrects;
    for i = 2:size(union_corrects,1)
        union_corrects(i,:) = max(union_corrects(i-1:i,:));
    end


    total_corrects = sum(union_corrects,2)';
    percentage_of_hits(i_file,:) = 100*total_corrects/total_frames;

    % Para calculo del resultado global sobre todos los archivos.
    total_global = total_global + total_frames;
    corrects_global = corrects_global + total_corrects;
end 

percentage_corrects_global = 100*corrects_global./total_global;
save(results_filename, 'percentage_of_hits','percentage_corrects_global')

percentage_corrects_global

%% Graficas

figure
x_tick = 1:num_db_files;
bar(x_tick, percentage_of_hits)
axis tight;
ax1 = gca;
set(ax1,'XTick',x_tick,'XTickLabel',names_cell)
title('Percentage of first candidate with correct f0 value.')
xlabel('Filename')
axis([x_tick(1)-0.5 x_tick(end)+0.5 0 100]);
colormap summer
