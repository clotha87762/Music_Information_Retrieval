%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function  [f0_medias chirp_rates t true_pitch_melody tsamples inds_present f0_value f0_time] = chirp_rate_estimation(fs, Nsamps, Lsamps, hop, labels_params)
%function  [f0_medias chirp_rates t true_pitch_melody tsamples inds_present f0_value f0_time] = chirp_rate_estimation(fs, Nsamps, Lsamps, hop, labels_params)

% time sample indexes of each signal frame
frame_time_samples = (1:hop:(Lsamps-Nsamps));
% sample time instants to interpolate lables (can be less dense!)
tsamples = (1:1:Lsamps)/fs;
% frame indexes
frame_indexes = (1:Nsamps);
% frame time instants
t = (frame_time_samples + Nsamps/2)/fs;
% number of frames
nframes = length(frame_time_samples);

[f0_value f0_time] = read_f0_labels(labels_params.database_path, labels_params.audio_file, labels_params.database);
% melody label interpolation in sample time instants
min_valid_f0 = 40; % Hz
true_pitch = interp1(f0_time, f0_value, tsamples);
valid_pitch = interp1(f0_time, f0_value > min_valid_f0, tsamples) == 1;
true_pitch = true_pitch .* valid_pitch;
inds_present = find(valid_pitch == 1);
true_pitch_melody = true_pitch(inds_present);

% median f0 estimated for each frame
f0_medias = zeros(nframes,1);
% chirp_rate estimated for each frame
chirp_rates = zeros(nframes,1);

% vector to compute weighted median (hanning window normalized)
Vm = hanning(Nsamps);
Vm = Vm / norm(Vm);
% vector to compute weighted slope (slope weighted by a hanning window)
Vpo = ((1:Nsamps)' - Nsamps/2)/fs;
Vp = Vpo .* (hanning(Nsamps).^2);
Vp = Vp - Vm * (Vp' * Vm); % in order Vp and Vm to be orthogonal
Vp = Vp / norm(Vp);

for ind = 1:nframes
    % sample indexes in the frame
    ind_samples = frame_time_samples(ind)+frame_indexes-1;
    % number of valid f0 values in the frame
    num_valid_pitch = sum(valid_pitch(ind_samples));
    % consider frame if at least Nsamps/4 have f0 values
    if (num_valid_pitch > Nsamps/4)
        % if whole frame has f0 values
        if (num_valid_pitch == Nsamps)
            Vml = Vm;
            Vpl = Vp;
        else
            % consider a masking function given by the valid pitches
            % this provides a correct estimation near boundaries
            Vml = Vm .* valid_pitch(ind_samples)';
            Vml = Vml / norm(Vml);
            Vpl = Vp .* valid_pitch(ind_samples)';
            Vpl = Vpl - Vml * (Vpl' * Vml);
            Vpl = Vpl / norm(Vpl);
        end
        
        % estimations
        f0_medias(ind) = true_pitch(ind_samples) * Vml / sum(Vml);
        chirp_rates(ind) = ((true_pitch(ind_samples)-f0_medias(ind)) / f0_medias(ind)) * Vpl / (Vpl'*Vpo);
        
%         figure(200)
%         %plot((true_pitch(ind_samples) / f0_medias(ind)) / (Vml'*Vmo), 'k')
%         plot((true_pitch(ind_samples) / f0_medias(ind)), 'k')
%         hold on
%         plot(1+chirp_rates(ind)*Vpo, 'r')
%         hold off
%         pause(0.01)
    end
end


end
