  function [spec t f f0gram f0s f0_hyps_indxs val_f0_hyps selected_warps] = stfcht(audio_file)
% function [spec t f f0gram f0s f0_hyps_indxs val_f0_hyps selected_warps] = stfcht(audio_file)
% 
% Function to compute a Short Time Fan Chirp Transform (STFChT) as well as
% a F0gram and pitch candidates for each frame. Parameters for the analysis 
% are described and must be set in the set_parameters.m file. 
% 

% set FChT parameters
set_parameters;

% add path to the c_code and m_code directories
addpath(['.' dir_slash 'c_code']);
addpath(['.' dir_slash 'm_code']);

% read audio file
[y fs] = audioread(audio_file);
% only one channel (left)
if (size(y,2) > 1); y = y(:,1); disp('WARNING: only left channel considered'); end

% audio filename to derive labels filename 
labels_params.audio_file = audio_file;

% design for efficient fcht computation
[warps cqt f0s accums] = design_fcht(fs, nfft, fmax, cqt_params, warp_params, f0_params, glogs_params);

% fcht computation
tic
[f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec selected_warps] = compute_fcht(y, fs, nfft, fmax, hop, cqt, warps, f0s, accums, f0_params, glogs_params, labels_params);
toc

ind_fmax = ceil(fmax/warps.fs_warp*nfft);
f = linspace(0,fmax,ind_fmax+1);

end
