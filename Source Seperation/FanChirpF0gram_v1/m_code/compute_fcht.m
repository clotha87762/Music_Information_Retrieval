 function [f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec selected_warps f0_medians chirp_rates tchr] = compute_fcht(y, fs, nfft, fmax, hop, cqt, warps, f0s, accums, f0_params, glogs_params, labels_params)
%function [f0gram t ind_maxs f0_hyps_indxs val_f0_hyps spec selected_warps f0_medians chirp_rates tchr] = compute_fcht(y, fs, nfft, fmax, hop, cqt, warps, f0s, accums, f0_params, glogs_params, labels_params)
%
% computes a spectrogram based on the FChT, a f0gram and predominant f0 hypotesis
%
%  in: 
%               y - audio signal
%              fs - sampling frequency
%            nfft - number of fft points
%            fmax - maximum frequency of interest (Hz)
%             hop - hop in samples
%
%             cqt - constant-Q transform parameters, struct with fields:
%                    .use_CQT        - whether to compute a constant-Q FChT 
%                    .Q              - quality factor Q for the constant-Q
%                    .db             - gain dB drop to define window with
%                    .softness       - smoothnes of the window width truncation
%                    .max_frac_pi    - determines maximum window length  
%                    .cqt.prepad     - pre padding value
%                    .cqt.pospad     - pos padding value
%                    .cqt.poles      - poles of the IIR LTV filter
%
%                    For a detailed description of this parameters see:
%                    "An efficient multi-resolution spectral transform for music analysis"
%                    P. Cancela, M. Rocamora, E. López, ISMIR 2009. Kobe, Japan.
%
%           warps - warpings parameters design
%                   data structure to do the interpolation using "c" code
%                   with fields:
%                    .pos1.int       - "c" index of previous sample
%                    .pos1.frac      - fractional value 
%                    .fs_orig        - sampling frequency after oversampling
%                    .fs_warp        - sampling frequency of warped signal
%                    .nsamps_torig   - number of samples in the oversampled signal frame
%                    .fact_over_samp - oversampling factor
%                    .chirp_rates    - chirp rate values
%
%             f0s - f0 grid values where to compute the f0gram
%
%          accums - spectrum array structs for efficient computation of the
%                   log spectrum accumulation, efficient attenuation of 
%                   mutliples and submultiples, and high pass filtering 
%
%                   .harmonic_accumulations: an array of structs with fields,
%                           .pos_int:       array of integer indexes for fast interpolation
%                           .pos_frac:      array of fractionary part for fast interpolation
%                           .weights:       array of weights to be aplied to addends
%                           .accum_length:  array of lengths of accumulations
%                           .global_weight: array of global weights applied to the accumulations
%
%                   .harmonic_attenuations: struct with the same fields as above
%
%                   .HP_logS: if high pass filtering is applied (1:yes 0:no)
%                   .HP_logS_components: values used in high pass filtering
%
%       f0_params - f0-gram parameters, struct with fields:
%                    .f0_params.f0min           - minimun fundamental frequency
%                    .f0_params.num_octs        - number of octaves of the f0-gram
%                    .f0_params.num_f0s_per_oct - number of f0s per octave 
%                    .f0_params.num_f0_hyps     - number of f0 hypotesis to extract
%                    .f0gram_type               - type of f0gram to return:
%                                                  'none' - no f0gram is returned (default)
%                                                'single' - single optimum direction f0gram
%                                              'compound' - best direction for each f0
%                    .f0_params.save_f0gram     -  whether to save the f0gram
%
%    glogs_params - gathered log spectrum parameters, struct with fields:
%                    .HP_logS             - whether to high-pass the GlogS
%                    .att_subharms        - whether to attenuate subharmonics
%                    .apply_GlogS_model   - whether to apply the normalization model
%                    .compute_GlogS_model - whether to learn the normalizartion model
%                    .median_poly_coefs   - model parameter value (median correction)
%                    .sigma_poly_coefs    - model parameter value (sigma  correction)
%
%   labels_params - parameters related to use of labels of fundamental
%                   frequency as the ones provided in the RWC or MIREX
%                   melody extraction test set, struc with fields:
%                   .chirp_rate_from_labels - use best chirp rate according to f0 labels
%                   .only_melody_frames     - process only frames labeled as melody
%                   .database               - name of the database ('mirex','rwc')
%                   .database_path          - path to the database files
%                   .audio_file             - audio filename (used to set labels filename)
%                   .fs_database            - sampling frequency of the database files
%                   .path_mat               - path to save mat files with the computations
%
% out: 
%
%          f0gram - F0gram, salience value for each f0 of the grid along audio frame instants
%               t - central time instants of frames
%        ind_maxs - warping index of optimal warping for each frame
%   f0_hyps_indxs - indexes of the f0 candidates in the grid for each frame
%     val_f0_hyps - F0gram amplitude value for each f0 candidate for each frame
%            spec - spectrum for optimal warping for each frame, i.e. STFChT 
%  selected_warps - indexes of selected warps for each f0 hypotesis for each frame
% 
% if chirp rates are estimated from f0 labels (labels_params.chirp_rate_from_labels > 0):
%      f0_medians - median f0 estimated for each frame using labels
%     chirp_rates - chirp_rate estimated from labels for each frame
%            tchr - frame time instants of estimations from labels
%

if nargout > 5; save_spec = 1; else save_spec = 0; end

% number of warpings
num_warps = size(warps.pos1_int,2);
% if chirp rate is derived from labels
if(labels_params.chirp_rate_from_labels > 0)
    Lsamps =  length(y);
    Nsamps = 2048;
    hopsamps = hop / warps.fact_over_samp;
    [f0_medians chirp_rates tchr] = chirp_rate_estimation(fs, Nsamps, Lsamps, hopsamps, labels_params);
    chirp_rates(isnan(chirp_rates)) = 0;
    if labels_params.chirp_rate_from_labels == 1
        num_warps = 1;
    end    
end

% number of samples of the warped signal frame
nsamps_twarp = size(warps.pos1_int,1);
% windows applied to warped signal frames
windows = repmat(hanning(nsamps_twarp),[1,num_warps]);

% frame indexes
frame_indexes = (1:warps.nsamps_torig);
% time sample indexes of each signal frame
frame_time_samples = (1:hop:warps.fact_over_samp*(length(y)-warps.nsamps_torig));
% number of frames
nframes = length(frame_time_samples);
% number of f0 values
nf0s = length(f0s);
% frequency bin index corresponding to fmax
ind_fmax = ceil(fmax/warps.fs_warp*nfft);
% minimun relative distance between f0 hypotesis peaks
min_peaks_dist = 1/100;

% optimum f0gram
if ~strcmp(f0_params.f0gram_type, 'none')
    f0gram = zeros(nf0s, nframes, 'single');
else
    f0gram = [];
end
% warping index of optimal warping
ind_maxs = zeros(nframes,1);
% indexes of f0 hypotesis for each frame
f0_hyps_indxs = zeros(f0_params.num_f0_hyps,nframes);
% values of f0 hypotesis for each frame
val_f0_hyps = zeros(f0_params.num_f0_hyps,nframes);
% indexes of selected warps for each f0 hypotesis for each frame
selected_warps = zeros(f0_params.num_f0_hyps,nframes);
% spectrum for optimal warping
if save_spec; spec = zeros((ind_fmax+1), nframes, 'single'); end

% whether to compute multiple f0 hypotesis
if f0_params.num_f0_hyps > 1; compute_mult_f0_hyps = 1; else compute_mult_f0_hyps = 1; end

% ====== GLOGS CORRECTION MODEL ======
if glogs_params.apply_GlogS_model
    % fundamental frequency indexes
    pos = 1:length(f0s);
    % model correction
    glogs_params.median_correction = polyval(glogs_params.median_poly_coefs,pos)';
    glogs_params.sigma_correction  = polyval(glogs_params.sigma_poly_coefs, pos)'*polyval(glogs_params.sigma_poly_coefs,pos(1));
end

% ====== PREPROCESSING ======
% anti-aliasing filter
ord = 20;
[B,A] = butter(ord,fmax/(fs/2));
y = filtfilt(B,A,y);

% oversampling 
y = resample(y,warps.fact_over_samp,1);

% ======  PROCESSING  ======
% if zero padding is done
zeropadding = nfft > nsamps_twarp;

if zeropadding
    x_mat = zeros(nfft,num_warps);
end

% frame time instants
t = (frame_time_samples + warps.nsamps_torig/2)/warps.fs_orig;

% interpolate f0 labels at frame time instants
if labels_params.chirp_rate_from_labels == 1
    chirp_rates_interp = interp1(tchr, chirp_rates, t);
end

% set frame indexes where to compute f0 estimates (labeled frames or all)
if labels_params.only_melody_frames == 1
    % load labels to ignore not melody frames
    [f0_value f0_time] = read_f0_labels(labels_params.database_path, labels_params.audio_file, labels_params.database);
    % melody label interpolation in frame time instants
    min_valid_f0 = 40; % Hz
    valid_pitch = interp1(f0_time, f0_value > min_valid_f0, t);
    index_of_frames = find(valid_pitch == 1);
else
    index_of_frames = 1:1:nframes;
end

% compute f0 preference function, if needed
if f0_params.prefer == 1
    f0_preference_weights = repmat(f0_preference_weighting(f0s, f0_params),[1,num_warps]);
end

for pos = index_of_frames
    
    % signal frame
    i =frame_time_samples(pos);
    x_frame = y(i+frame_indexes,1);
    
    if labels_params.chirp_rate_from_labels == 1
        [val index_warps] = min(abs(chirp_rates_interp(pos)-warps.chirp_rates));
    else
        index_warps = 1:num_warps;    
    end
    
    % warped frames 
    x_warpings = quick_interp1_int_frac_c(x_frame, warps.pos1_int(:,index_warps), warps.pos1_frac(:,index_warps));
    
    % fft of warped frames (FChT)
    if cqt.use_CQT % constant Q FChT
        fcht = zeros(nsamps_twarp/2,num_warps);
        transformada_shift = fft(x_warpings);
        for j = 1:num_warps
           fcht(:,j) = iir_cqt_c(transformada_shift(:,j),cqt.poles,cqt.prepad,cqt.pospad);
        end  
        if zeropadding
            x_mat(1:nsamps_twarp,:) = fftshift(ifft([fcht; zeros(1,num_warps) ;conj(flipud(fcht(2:end,:)))]),1);
            fcht = fft(x_mat);       
        end
    else
        if zeropadding
            x_mat(1:nsamps_twarp,:) = x_warpings.*windows;
            fcht = fft(x_mat);
        else
            fcht = fft(x_warpings.*windows);
        end
    end    
    
    % consider spectrums only up to fmax
    fcht = fcht(1:(ind_fmax+1),:);
    
    % log spectrum
    logS = log10_1_alpha_abs_c(fcht,10); % logS = log10(1+10*abs(fcht));
      
    % harmonic accumulation to compute pitch salience
    f0gram_frame = harmonic_accumulation(num_warps, logS, accums.harmonic_accumulations, accums.harmonic_attenuations, nf0s, f0_params.num_f0s_per_oct, accums.HP_logS_components, glogs_params);

    % apply f0 preference weighting, if needed
    if f0_params.prefer == 1
        f0gram_frame = f0gram_frame .* f0_preference_weights;
    end
    
    % predominant pitch estimation
    [fool indx] = max(max(f0gram_frame));    
    ind_maxs(pos) = indx;
    
    % save optimum spectrum (if needed) and f0gram
    if save_spec; spec(:,pos) = fcht(:,indx); end
    
    % whether to return a single optimum direction f0gram
    if strcmp(f0_params.f0gram_type, 'single')    
        f0gram(:,pos) = f0gram_frame(:,indx);
    end
    % whether to compute compound f0gram and return it   
    if (strcmp(f0_params.f0gram_type, 'compound') || compute_mult_f0_hyps)   
        % optimum warping for each f0
        [f0gram_maxs chirp_rate_indexes] = max(f0gram_frame,[],2);
        if strcmp(f0_params.f0gram_type, 'compound')
            f0gram(:,pos) = f0gram_maxs;
        end
    end    
        
    % compute multiple f0 hypotesis
    if compute_mult_f0_hyps
        % find local maxima of compound f0gram
        ind_maxs = [];
        for j = 2:nf0s-1;
           if ((f0gram_maxs(j) > f0gram_maxs(j-1)) && (f0gram_maxs(j) > f0gram_maxs(j+1)))
               ind_maxs = [ind_maxs; j];
           end
        end
        % sort local maxima
        [vals,inds] = sort(f0gram_maxs(ind_maxs),'descend'); 
        % save biggest local maxima
        num_max_inds = length(inds);
        if num_max_inds > 0
            num_inds = min(f0_params.num_f0_hyps,num_max_inds);
            % add f0 candidates considering they are 1% appart
            f0_hyps_indxs(1,pos) = ind_maxs(inds(1));  
            val_f0_hyps(1,pos) = f0gram_maxs(f0_hyps_indxs(1,pos)); 
            if labels_params.chirp_rate_from_labels == 1
                selected_warps(1,pos) = warps.chirp_rates(index_warps);
            else
                selected_warps(1,pos) = warps.chirp_rates(chirp_rate_indexes(ind_maxs(inds(1))));
            end
            k = 2;
            for m = 2:num_inds
                fin = 0;
                while ((k <= num_max_inds) && ~fin) 
                    fin = (min(abs(f0s(f0_hyps_indxs(1:m-1,pos)) - f0s(ind_maxs(inds(k))))./f0s(f0_hyps_indxs(1:m-1,pos))) > min_peaks_dist);
                    if ~fin
                        k = k + 1;
                    end
                end
                if fin
                    f0_hyps_indxs(m,pos) = ind_maxs(inds(k));
                    val_f0_hyps(m,pos) = f0gram_maxs(f0_hyps_indxs(m,pos));
                    if labels_params.chirp_rate_from_labels == 1
                        selected_warps(m,pos) = warps.chirp_rates(index_warps);
                    else
                        selected_warps(m,pos) = warps.chirp_rates(chirp_rate_indexes(ind_maxs(inds(k))));
                    end
                    k = k + 1;
                end
            end
        end
    else
        [val f0_hyps_indxs(pos)] = max(f0gram_frame(:,indx));
    end

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accum = accumulate_vectors(x,accum_interpolations,num_f0s_per_oct)
% function accum = accumulate_vectors(x,accum_interpolations)
%
% accumulation of the signal x in position given in accum_interpolations
% 
% this is done to compute the salience of diferent f0 candidates, by
% accumulation of the spectrum at the position of harmonics
% each accumulation is then weighted by a global weigth that in this case
% corresponds to the number of accumulated points
%

% interpolation of values at postions given by accum_interpolations
interps = quick_interp1_int_frac_c(x,accum_interpolations.pos_int,accum_interpolations.pos_frac);

%accum = accum_arrays(interps,accum_interpolations.accum_length);
accum = recursive_octave_accum(accum_arrays(interps,accum_interpolations.accum_length)',num_f0s_per_oct)';

% weighting of accumulated values by the number of accumulated points
accum = accum.*accum_interpolations.global_weight;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f0gram_frame = harmonic_accumulation(num_warps, logS, accum_interpolations, subharmonic_interpolations, num_f0s, num_f0s_per_oct, HP_logS_components, glogs_params)
%function f0gram_frame = harmonic_accumulation(num_warps, logS, accum_interpolations, subharmonic_interpolations, num_f0s, num_f0s_per_oct, HP_logS_components, HP_logS, att_subharms, compute_GlogS_model, apply_GlogS_model, median_correction, sigma_correction)

% high pass filter and mean substraction
if(~isempty(HP_logS_components))
    if(glogs_params.HP_logS == 1)
        logS = logS - HP_logS_components*(logS'*HP_logS_components)';
    end
    for j = 1:num_warps
        if(glogs_params.HP_logS == 2)
            logS(:,j) = quick_filtfilt_1_v2010(B,A,logS(:,j));
        end
        logS(:,j) = (logS(:,j)-mean(logS(:,j)));
    end
end

% global variables to learn the GlogS model
if glogs_params.compute_GlogS_model
    global mean_glogs
    global mean_glogs_post
    global mean_glogs2
    global mean_glogs_post2
    global cant_N
end

f0gram_frame = zeros(num_f0s, num_warps);

for j = 1:num_warps
        
       % accumulation of log spectrum
       glogs = accumulate_vectors(logS(:,j),accum_interpolations, num_f0s_per_oct);        

       % post-processing of glogs for subharmonic attenuation
       f0gram_frame(:,j) = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, glogs_params.att_subharms, num_f0s_per_oct);              
       
       % apply model correction 
       if glogs_params.apply_GlogS_model
           %f0gram_frame(:,j) = (f0gram_frame(:,j)-polyval(P,pos)')./(polyval(Ps,pos)'*polyval(Ps,pos(1)));
           f0gram_frame(:,j) = (f0gram_frame(:,j)-glogs_params.median_correction)./glogs_params.sigma_correction;
       end
       
       % model learning
       if glogs_params.compute_GlogS_model       
           if(isempty(mean_glogs))
               mean_glogs = glogs;
               mean_glogs_post = f0gram_frame(:,j);
               mean_glogs2 = glogs.*glogs;
               mean_glogs_post2 = f0gram_frame(:,j).*f0gram_frame(:,j);
               cant_N = 1;
           else
               cant_N = cant_N + 1;
               mean_glogs = (mean_glogs*(cant_N-1) + glogs)/cant_N;
               mean_glogs_post = (mean_glogs_post*(cant_N-1) + f0gram_frame(:,j))/cant_N;
               mean_glogs2 = (mean_glogs2*(cant_N-1) + glogs.*glogs)/cant_N;
               mean_glogs_post2 = (mean_glogs_post2*(cant_N-1) + f0gram_frame(:,j).*f0gram_frame(:,j))/cant_N;
           end
       end

end
    

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glogs_pospro = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, att_subharms, num_f0s_per_oct)
% function glogs_pospro = post_process_glogs(glogs, subharmonic_interpolations, num_f0s, att_subharms, num_f0s_per_oct)
%
% postprocessing of the gathered log spectrum by attenuation of each value 
% by the maximum of the values at the subharmonics
% and subharmonics attenuations by the scaled value at the next octave

if ~att_subharms 
    glogs_pospro = glogs(end-num_f0s+1:end) - max_arrays(quick_interp1_int_frac_c(glogs,subharmonic_interpolations.pos_int,subharmonic_interpolations.pos_frac), subharmonic_interpolations.accum_length);
else
    glogs_pospro = glogs((num_f0s_per_oct+1):end) - max_arrays(quick_interp1_int_frac_c(glogs,subharmonic_interpolations.pos_int,subharmonic_interpolations.pos_frac), subharmonic_interpolations.accum_length);
    %glogs_pospro = glogs_pospro(1:num_f0s) - 1/2 * glogs_pospro(num_f0s_per_oct+1:end);
    glogs_pospro = glogs_pospro(1:num_f0s) - 1/3 * glogs_pospro(num_f0s_per_oct+1:end);
end
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function f0_preference_weights = f0_preference_weighting(f0s, f0_params)
% models f0 preference, for instance for melody detection, as a gaussian
% whose parameters are derived for example from labeled data

MIDI_grid = 69 + 12 * log2(f0s/440);

bell = 1/sqrt(2*pi*f0_params.prefer_stdev^2).*exp(-(MIDI_grid-f0_params.prefer_mean).^2/(2*f0_params.prefer_stdev^2));

% kind of simple smoothing for non observed frequencies
smoothing_offset = 0.01;
f0_preference_weights = ( smoothing_offset + bell') / ( smoothing_offset + 1);

end
