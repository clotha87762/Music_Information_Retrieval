% test stfcht
% 
% Example code to show how to compute a spectrogram based on the FChT as
% well as a F0gram, and f0 candidates. Parameters for the analysis are 
% described and must be set in the set_parameters.m file. 
% 
clear all
close all

%% AUDIO FILE
audio_file = 'c5.wav';

%% STFChT
[spec t f f0gram f0s f0_hyps_indxs val_f0_hyps selected_warps] = stfcht(audio_file);

%% PLOTS

A4 = 440; % reference A4 frequency for pitch grid in plots
labels_on = 0; % whether to use f0 labels in plots
chirps_from_labels = 0; % whether to plot chirp rates estimated from f0 labels

num_semitones_down = floor(12*log2(A4/f0s(1)));
num_semitones_up   = floor(12*log2(f0s(end)/A4));
value_ticks = A4*2.^((-num_semitones_down:num_semitones_up)/12);
pos_ticks = interp1(f0s, 1:length(f0s), value_ticks);
f_ylabels = num2str(value_ticks', '%4.2f');
black = abs(1-gray);

notes_names = {'A-'; 'Bb'; 'B-'; 'C-'; 'C#'; 'D-'; 'Eb'; 'E-';  'F-'; 'F#'; 'G-'; 'Ab'} ;
semitones = (-num_semitones_down:num_semitones_up)';
octaves = num2str(floor((semitones+9)/12) + 4,'%1d');
semitones_notes_names = notes_names(mod(semitones,12)+1);
ylabels = semitones_notes_names;
for i=1:length(semitones)
    ylabels{i} = strcat(ylabels{i}, octaves(i));
    ylabels{i} = strcat(ylabels{i}, '-');
    ylabels{i} = strcat(ylabels{i}, f_ylabels(i,:));
end

% plot of Short Time Fan Chirp Transform
imagesc(t,f,20*log10(double(abs(spec))+0.1))
colormap(black);set(gca,'YDir','normal');
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('STFChT: spectrogram for a single optimum warping per frame')

% plot of F0gram
figure
hh = imagesc(t,1:length(f0s),double(f0gram));
colormap(black);set(gca,'YDir','normal');
title('F0gram for the best direction for each f0')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(get(hh,'Parent'),'YTick',pos_ticks)
set(get(hh,'Parent'),'YTickLabel',ylabels)
ylim([1 length(f0s)])
grid


% load labels to display and compute chirp rates 
if labels_on == 1
    % set FChT parameters
    set_parameters;
    % audio filename to derive labels filename 
    labels_params.audio_file = audio_file;
    % read lables
    [f0_value f0_time] = read_f0_labels(labels_params.database_path, labels_params.audio_file, labels_params.database);
    f0_value_ind = interp1(f0s,1:length(f0s),f0_value);
    
    % whether to estimate chirp rates from labels
    if chirps_from_labels == 1
        % read audio file
        [y fs] = wavread(audio_file);
        Lsamps =  size(y,1);
        Nsamps = 2048;
        hopsamps = hop / warp_params.fact_over_samp;
        [f0_medians chirp_rates tchr] = chirp_rate_estimation(fs, Nsamps, Lsamps, hopsamps, labels_params);
        chirp_rates(isnan(chirp_rates)) = 0;
        if labels_params.chirp_rate_from_labels == 1; num_warps = 1; end    
    end
end

% same plot with f0 candidates
hop = 3;
f0_hyps = f0_hyps_indxs(:,1:hop:end);
t_f0_hyps = t(1:hop:end);

figure
hh = imagesc(t,1:length(f0s),double(f0gram));
colormap(black);set(gca,'YDir','normal');
title('F0gram and f0 candidates')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(get(hh,'Parent'),'YTick',pos_ticks)
set(get(hh,'Parent'),'YTickLabel',ylabels)
hold on
plot(t_f0_hyps,f0_hyps(1,:),'ok','MarkerSize',4)
plot(t_f0_hyps,f0_hyps(2,:),'+k','MarkerSize',4,'LineWidth',0.1)
plot(t_f0_hyps,f0_hyps(3,:),'xk','MarkerSize',4,'LineWidth',0.1)
legend('First','Second','Third');
if labels_on == 1
    plot(f0_time,1.03*f0_value_ind,'k','LineWidth',0.05)
    plot(f0_time,0.97*f0_value_ind,'k','LineWidth',0.05)
    title('F0gram and f0 candidates (and 3% band centered at f0 label)')
end
hold off
ylim([1 length(f0s)])
grid


% selected chirp rates 
figure;
plot(t,selected_warps(1,:),'ok', 'MarkerSize',4);
if chirps_from_labels == 1
    hold on;
    plot(tchr,chirp_rates,'-xk','Color',[0.6,0.6,0.6], 'MarkerSize',4);
    hold off
    legend('selected chirp rates', 'chirp rates from labels')    
end
xlabel('Time (s)')
ylabel('\alpha')
title('Chirp rate estimation for best candidate')
axis tight
grid
