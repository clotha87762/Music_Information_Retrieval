%% main_visualize_individual_results.m
%
% Program to visualize the pitch candidates estimation and the f0 labels.
% The F0gram and the pitch candidates extraction must be previously done.

close all;
clear all;

% set FChT parameters
set_parameters;

% audio filename
labels_params.audio_file = 'daisy4';

% mat files of the analysis results
f0gram_filename = [labels_params.path_mat labels_params.audio_file '_f0gram.mat'];
pitch_filename = [labels_params.path_mat labels_params.audio_file '_pitch.mat'];

% load results files
if f0_params.save_f0gram; 
    load(f0gram_filename); 
end
load(pitch_filename);

% load labels
[f0_value f0_time] = read_f0_labels(labels_params.database_path, labels_params.audio_file, labels_params.database);

%% PLOTS

black = abs(1-gray);
colors = 'bcgkym';
legends = ['1st';'2nd';'3rd';'4th';'5th';'6th'];    
num_f0_hyps = min(size(f0_hyps_indxs,1), length(colors));


if f0_params.save_f0gram
    figure
    imagesc(t,1:length(f0s),f0gram)
    colormap(black);set(gca,'YDir','normal');
    title('F0gram')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    hold on
    for i = 1:num_f0_hyps
        plot(t,f0_hyps_indxs(i,:),['.' colors(i)])
    end    
    hold off
    legend(legends(1:num_f0_hyps,:))
end

f0_hyps_indxs(f0_hyps_indxs == 0) = 1;
figure;
plot(f0_time,f0_value,'.r') 
hold on
for i = 1:num_f0_hyps
    plot(t,f0s(f0_hyps_indxs(i,:)),['.' colors(i)])
end    
hold off
axis tight
legend(['f_0'; legends(1:num_f0_hyps,:)])
title('f0 hypotesis and label')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
