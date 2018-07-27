%% main_learn_normalization_model
%
% Program to process a database of audio files computing the FChT and the F0gram
% for each file. The mean and the variance of the Gathered Log Spectrum (GlogS)
% are collected for each f0 for every frame of the music collection (e.g. RWC). 
% This is used to normalize the GlogS to obtain zero median and unit variance. 
% The normalization is done by means of a polynomial approximation to the
% collected statistics for the median and the variance. The polynomials
% evaluated at each f0 are the model used to obtain the normalization. 

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

% make sure NOT to apply the GlogS model correction
glogs_params.apply_GlogS_model = 0;
% make sure to learn the GlogS model
glogs_params.compute_GlogS_model = 1; % 0 = no, 1 = global, 2 = individual
% filename to save glogs statistics to build the model
glogs_params.model_filename = ['glogs_model_' labels_params.database '.mat'];

% design for efficient fcht computation
[warps cqt f0s accums] = design_fcht(labels_params.fs_database, nfft, fmax, cqt_params, warp_params, f0_params, glogs_params);

% variables for learning the GlogS normalization model
global mean_glogs
global mean_glogs_post
global mean_glogs2
global mean_glogs_post2
global cant_N


%% ============= DATABASE FILES =============

% struct with the name and complete path to each audio file in the database 
db_files = files_in_path(labels_params.database_path, 'wav');
% number of files in the database
num_db_files = length(db_files);


%% ====== PROCESSING ======

tic
for i_file = 1 : num_db_files

    if(glogs_params.compute_GlogS_model == 2)
         mean_glogs = 0;
         mean_glogs_post = 0;
         mean_glogs2 = 0;
         mean_glogs_post2 = 0;
         cant_N = 0;
    end
    
    % audio filename without path
    labels_params.audio_file = db_files(i_file).clean_file_name;

    % mat files for the result of the analysis
    f0gram_filename = [labels_params.path_mat labels_params.audio_file '_f0gram.mat'];
    pitch_filename = [labels_params.path_mat labels_params.audio_file '_pitch.mat'];

    % if pitch detection has NOT been done previously (file does not exist)
    if ~(exist(pitch_filename,'file') == 2)

        % audio file reading
        disp(['Processing file ' labels_params.audio_file ' ...'])
        [y, fs] = wavread(db_files(i_file).name);
        if (labels_params.fs_database ~= fs); disp('WARNING: audio fs different from design fs'); end;
        % only one channel (left)
        if (size(y,2) > 1); y = y(:,1); disp('WARNING: only left channel considered'); end
        
        % fcht computation
        [f0gram t ind_maxs f0_hyps_indxs val_f0_hyps] = compute_fcht(y, fs, nfft, fmax, hop, cqt, warps, f0s, accums, f0_params, glogs_params, labels_params);
        
        % save results
        save(pitch_filename, 'f0s','f0_hyps_indxs','val_f0_hyps','t');
        if f0_params.save_f0gram; save(f0gram_filename, 'f0gram','t','f0s'); end
                
        % compute model parameters
        sigma = 2*sqrt(mean_glogs2-mean_glogs.^2);
        sigma_post = 2*sqrt(mean_glogs_post2-mean_glogs_post.^2);
        % save the GlogS model
        if (glogs_params.compute_GlogS_model == 1)
            save(glogs_params.model_filename, 'mean_glogs', 'mean_glogs2', 'mean_glogs_post', 'mean_glogs_post2', 'sigma', 'sigma_post', 'cant_N', 'f0s');
        elseif (glogs_params.compute_GlogS_model == 2)
            glogs_individual_model_filename = ['glogs_model.' labels_params.audio_file '.mat'];    
            save(glogs_individual_model_filename, 'mean_glogs', 'mean_glogs2', 'mean_glogs_post', 'mean_glogs_post2', 'sigma', 'sigma_post', 'cant_N', 'f0s');
        end
        
        % plots of the GlogS_model
        figure(400)
        plot(mean_glogs)
        hold on
        plot(mean_glogs + sigma,'r')
        plot(mean_glogs - sigma,'r')
        hold off
        figure(401)
        plot(mean_glogs_post)
        hold on
        plot(mean_glogs_post + sigma_post,'r')
        plot(mean_glogs_post - sigma_post,'r')
        hold off
        pause(0.001)

    else 
        disp(['Skipping file ' labels_params.audio_file ' (remove .mat to force computation)'])
    end
    
end
toc 

%% estimate GlogS polinomial model aproximation 
if glogs_params.compute_GlogS_model == 1
    load(glogs_params.model_filename);
    
    pos = 1:length(f0s);
    format long
    disp('Model polynomial coefficients: ')
    P = polyfit(pos,mean_glogs_post',2)
    Ps = polyfit(pos,sigma_post',2)
        
    figure    
    plot(pos,mean_glogs_post,'b',pos,polyval(P,pos),'r')
    figure
    plot(pos,sigma_post,'b',pos,polyval(Ps,pos),'r')
    
    save(glogs_params.model_filename, 'P', 'Ps', 'pos','mean_glogs', 'mean_glogs2', 'mean_glogs_post', 'mean_glogs_post2', 'sigma', 'sigma_post', 'cant_N', 'f0s')
end
