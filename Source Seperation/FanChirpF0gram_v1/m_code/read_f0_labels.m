%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f0_value f0_time] = read_f0_labels(path, audio_file, database)
% function [f0_value f0_time] = read_f0_labels(path, audio_file, database)
%
% read f0 labels for the corresponding audio file given in a text file
%
%  in: 
%            path - full path to the audio file
%      audio_file - audio file name (without extension)
%        database - process the labels file according to a certain database
%                   'mirex': Mirex 2004 and 2005 Melody Contest files,
%                     'rwc': RWC Music Database melody files
%
% out: 
%        f0_value - f0 values for the time instants given in f0_time
%         f0_time - time instants of the labels 
%

if ~isempty(strfind(audio_file, '.wav'))
    ind_point = find(audio_file == '.', 1, 'last');
    audio_file = audio_file(1:ind_point-1);
end

if strcmp(database, 'mirex')
    
    reference_filename = [path audio_file 'REF.txt'];
    f0_ref = dlmread(reference_filename);
    f0_time = f0_ref(:,1);
    f0_value = f0_ref(:,2);
    
elseif strcmp(database, 'rwc')
    if(isunix); dir_slash = '/'; else dir_slash = '\'; end    
    reference_filename = [path dir_slash 'MELODY' dir_slash audio_file '.MELODY.TXT'];
    [num_frames fool1 fool2 pitch_value fool3] = textread(reference_filename, '%d %d %s %f %f');
    % RWC parameters
    hop_time_lbs = 0.01; N_lbs = 2048; fs = 44100; 
    nframes_lbs = num_frames(end)+1; 
    f0_value = zeros(nframes_lbs,1);
    f0_time  = ((1:nframes_lbs)-1) * hop_time_lbs + (N_lbs/2) / fs;
    f0_value(num_frames+1) = pitch_value; 
    
end


end