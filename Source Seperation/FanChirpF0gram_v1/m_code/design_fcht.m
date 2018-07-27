  function [warps cqt f0s accums] = design_fcht(fs, nfft, fmax, cqt_params, warp_params, f0_params, glogs_params)
% function [warps cqt f0s accums] = design_fcht(fs, nfft, fmax, cqt_params, warp_params, f0_params, glogs_params)
%
% design for efficient fcht computation
%
%  in:
%              fs - sampling frequency
%            nfft - number of fft points
%            fmax - maximum frequency of interest (Hz)
%
%      cqt_params - constant-Q transform parameters, struct with fields:
%                    .use_CQT        - whether to compute a constant-Q FChT 
%                    .Q              - quality factor Q for the constant-Q
%                    .db             - gain dB drop to define window with
%                    .softness       - smoothnes of the window width truncation
%                    .max_frac_pi    - determines maximum window length  
%
%                    For a detailed description of this parameters see:
%                    "An efficient multi-resolution spectral transform for music analysis"
%                    P. Cancela, M. Rocamora, E. López, ISMIR 2009. Kobe, Japan.
%
%     warp_params - warpings parameters
%                   struct with fields (may vary according to .warping_type)
%                    .warping_type   - warping type (linear/not-linear)
%                    .fact_over_samp - oversampling factor
%                    .num_warps      - number of warpings
%                    .nsamps_twarp   - number of samples in warped signal frame
%                    .alpha_max      - maximum normalized frequency deviation (linear)
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
% out: 
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
%            cqt - constant-Q transform parameters, struct with fields:
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

% ============= CONSTANT Q DESIGN =============

% number of points of pre and post padding used to set initial conditions
cqt.use_CQT = cqt_params.use_CQT;
nfft_cqt = warp_params.nsamps_twarp;
if cqt_params.use_CQT
    % pre and pos padding values
    cqt.prepad = 10;
    cqt.pospad = 20;
    % frequency bins grid (linear in this case) - pre and pos padding is added
    slope = pi/nfft_cqt;
    thetas = abs(slope*linspace(-cqt.prepad,nfft_cqt/2+cqt.pospad-1,nfft_cqt/2+cqt.prepad+cqt.pospad));
    thetas(cqt.prepad+1) = eps; % zero digital frequency
    % poles of the IIR LTV Q FFT transform for the parameters above
    cqt.poles = calculate_poles(thetas,cqt_params,nfft_cqt);
end

% =============  WARPING DESIGN   =============
% sampling frequency after oversampling
fs_orig = warp_params.fact_over_samp * fs;
% number of samples of the original signal frame
nsamps_torig = 4 * warp_params.fact_over_samp * warp_params.nsamps_twarp;
% time instants of the original signal frame
t_orig = ((1:nsamps_torig)-nsamps_torig/2)/fs_orig;

% warping design according to different warping types
if strcmp(warp_params.warping_type, 'linear')
    
    % linear chirps warping definition as relative frequency deviation 
    [freq_relative chirp_rates] = define_warps_linear_chirps(warp_params.alpha_max,warp_params.num_warps,t_orig);
    
    % maximum relative frequency deviation
    freq_relative_max = max(max(freq_relative));

    % sampling frequency of warped signal to be free of aliasing up to fmax
    fs_warp = 2 * fmax * freq_relative_max;

    % time instants of the warped signal frame
    t_warp = ((1:warp_params.nsamps_twarp)-warp_params.nsamps_twarp/2)/fs_warp;

    % design of warpings for efficient interpolation
    warps = design_warps(freq_relative, t_orig, t_warp, fs_orig, cqt_params.use_CQT);

    % return sampling frequency of warped signal 
    warps.fs_warp = fs_warp;

    % return sampling frequency after oversampling 
    warps.fs_orig = fs_orig;
    
    % return number of samples of the original signal frame
    warps.nsamps_torig = nsamps_torig;
    
    % return oversampling factor 
    warps.fact_over_samp = warp_params.fact_over_samp;

    % return chirp rates
    warps.chirp_rates = chirp_rates;
        
elseif strcmp(warp_params.warping_type, 'not-linear')
   disp('Error: Warping type not implemented') 
else
   disp('Error: Warping type not known') 
end

% ======== GATHERED LOG SPECTRUM DESIGN =======

% maximun fundamental frequency
if glogs_params.att_subharms
    f0max = f0_params.f0min * 2^(f0_params.num_octs+1);
else
    f0max = f0_params.f0min * 2^f0_params.num_octs;
end
% minimun fundamental frequency for harmonic attenuation
f0min_subh = f0_params.f0min/2;
% number of f0 values considering extra octaves for harmonic and subharmonic attenuation
if glogs_params.att_subharms
    % number of f0 values considering two extra octaves
    num_f0s = (f0_params.num_octs+2) * f0_params.num_f0s_per_oct + 1;
else
    % number of f0 values considering an extra octave
    num_f0s = (f0_params.num_octs+1) * f0_params.num_f0s_per_oct + 1;
end
% f0 values considering an extra octave for subharmonic attenuation
f0s = (logspace(log10(f0min_subh),log10(f0max),num_f0s));

% spectral bin positions and weights of harmonics for each f0
harmonic_accumulations_design = design_harmonic_accumulations(f0s,fs_warp,nfft,fmax,f0_params.num_f0s_per_oct);
% an array of structs to do the acummulation of the spectrum efficiently
harmonic_accumulations = design_weighted_accum_interpolations(harmonic_accumulations_design);

% extra octave is removed
if glogs_params.att_subharms
    f0s = f0s((f0_params.num_f0s_per_oct+1):end-f0_params.num_f0s_per_oct);
else
    f0s = f0s((f0_params.num_f0s_per_oct+1):end);
end

% harmonic attenuation design
harmonic_positions = design_harmonic_attenuation(f0_params.num_f0s_per_oct, num_f0s);
% an array of structs to do the attenuation of harmonics efficiently
harmonic_attenuations = design_weighted_accum_interpolations(harmonic_positions);

% return spectrum accumulation design
accums.harmonic_accumulations = harmonic_accumulations;
accums.harmonic_attenuations = harmonic_attenuations;

ind_fmax = ceil(fmax/warps.fs_warp*nfft);

accums.HP_logS_components = [];

accums.HP_logS = glogs_params.HP_logS;
if glogs_params.HP_logS  
    if(glogs_params.HP_logS ==1)        
        accums.HP_logS_components = linspace(-1,1,(ind_fmax+1))';
        accums.HP_logS_components = accums.HP_logS_components/norm(accums.HP_logS_components);
    end
    if(glogs_params.HP_logS ==2)
        [B,A] = cheby1(1,0.1,0.01,'high');
        accums.HP_logS_components.A = A;
        accums.HP_logS_components.B = B;
        accums.HP_logS_components.weights = linspace(1,0.1,(ind_fmax+1))';        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accumulations = design_harmonic_accumulations(f0s,fs,nfft,max_freq,n_octave)
% function accumulations = design_harmonic_accumulations(f0s,fs,nfft,max_freq)
%
% computes frequency positions of harmonics for each fundamental frequency
% which are used in the computation of the salience of different pitches
% these are stored in a structure along with a weight value for each
% position that can be used to attenuate or stress some parts of the spectrum
%
%  in:
%                f0s - fundamental frequency values
%                 fs - sampling frequency
%               nfft - number of fft points
%           max_freq - maximum frequency of interest
%
% out:
%      accumulations - an array of structs with fields: 
%         * .positions:    frequency position in bins of each harmonic
%         * .weigths:      weight value for each harmonic 
%         * .g_weigth:     a global weight
%
 
for i = 1:length(f0s)
    % frequency bin number of each harmonic
    if(nargin<5)
    partials_pos = (f0s(i):f0s(i):max_freq)/(fs)*(nfft)+1;
    else
    if(i>length(f0s)-n_octave)
        partials_pos = (f0s(i):f0s(i):max_freq)/(fs)*(nfft)+1;
    else
        partials_pos = (f0s(i):2*f0s(i):max_freq)/(fs)*(nfft)+1;
    end 
    end
    accumulations(i).positions = partials_pos;
    % weight value for each harmonic 
    accumulations(i).weights = ones(size(partials_pos))/length(partials_pos);
    accumulations(i).g_weight = 1/length((f0s(i):f0s(i):max_freq));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function accumulations = design_harmonic_attenuation(numf0s_oct, num_f0s)
% function accumulations = design_harmonic_attenuation(numf0s_oct, num_f0s)
%
% computes frequency positions of subharmonics for each fundamental frequency
% which are used in the computation of the harmonic attenuation of the
% salience function 
% these are stored in a structure along with a weight value for each
% position that can be used to attenuate or stress some parts
% 
%
%  in:
%         numf0s_oct - number of f0 values per octave
%             numf0s - total number of f0 values
%
% out:
%      accumulations - an array of structs with fields: 
%         * .positions:    frequency position in bins of each harmonic
%         * .weigths:      weight value for each harmonic 
%         * .g_weigth:     a global weight
%

% maximum number of subharmonics
ind_max_sh = floor(2^(num_f0s/(numf0s_oct)));
indx = 2:ind_max_sh;

for i = numf0s_oct+1:num_f0s
    
    pos_sh = i - (numf0s_oct * log2(indx));
    accumulations(i-numf0s_oct).positions = pos_sh(pos_sh >= 1);
    accumulations(i-numf0s_oct).weights = ones(size(accumulations(i-numf0s_oct).positions));
    accumulations(i-numf0s_oct).g_weight =1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function warps = design_warps(freq_relative, t_orig, t_warp, fs_orig, zero_phase_shift)
% function warps = design_warps(freq_relative, t_orig, t_warp, fs_orig)
% 
% design of a data structure to do the warping efficiently using "c" code
%
% the warping is done by interpolating the original signal in time instants
% given by the desired frequency deviation, to do this, the interpolation 
% instants are stored in a structure as an integer index and a fractional value
%
%  in:
%      freq_relative : relative frequency deviations
%             t_orig : original time vector
%             t_warp : warped   time vector
%            fs_orig : original sampling frequency
%
% out:
%              warps : data structure to do the interpolation using "c" code
%                      pos1.int  - "c" index of previous sample
%                      pos1.frac - fractional value 
%
% hypothesis: sampling frequency at the central point equals the original
%

% number of warpings
num_warps = size(freq_relative,2);

% interpolation points
pos1 = zeros(num_warps, length(t_warp));

% for each warping interpolation points are computed
for i = 1:num_warps
    if(min(freq_relative(:,i))<0)
        disp('Error: relative warping frequencies must be non-negative.');
        disp('       Setting relative frequencies to 1 for this warping.')
        freq_relative(:,i) =1;
    end
    % integration of relative frequency to obtain phase values
    phi = cumtrapz(t_orig,freq_relative(:,i)');
    % centering of phase values to force original frequency in the middle
    phi = phi - phi(end/2);
    % interpolation of phase values to obtain warped positions
    pos1(i,:) = interp1(phi,t_orig,t_warp)*fs_orig + length(t_orig)/2;
    if(zero_phase_shift)
        pos1(i,:) = fftshift(pos1(i,:));
    end
        
end

% previous sample index
pos1_int = uint32(floor(pos1))';
% integer corresponding to previous sample index in "c"
warps.pos1_int = (pos1_int - uint32(1));
% fractional value that defines the warped position 
warps.pos1_frac = (double(pos1)' - double(pos1_int));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  acum_interpolations = design_weighted_accum_interpolations(accumulations)
% function  acum_interpolations = design_weighted_accum_interpolations(accumulations)
%
% design of a data structure to do the accumulation of spectrum efficiently using "c" code
%
%  in:
%           accumulations - an array of structs with fields: 
%         * .positions:    frequency position in bins of each harmonic
%         * .weigths:      weight value for each harmonic 
%
% out:
%     acum_interpolations - an array of structs with fields: 
%         * .pos_int:       array of integer indexes for fast interpolation
%         * .pos_frac:      array of fractionary part for fast interpolation
%         * .weights:       array of weights to be aplied to addends
%         * .accum_length:  array of lengths of accumulations
%         * .global_weight: array of global weights applied to the accumulations

% memory allocation for the array of structs for effiency
cant_addends = length([accumulations(:).positions]);
acum_interpolations.pos_int = uint32(zeros(cant_addends,1));
acum_interpolations.pos_frac = zeros(cant_addends,1);
acum_interpolations.weights = zeros(cant_addends,1);
acum_interpolations.accum_length = length(accumulations);
acum_interpolations.global_weight = length(accumulations);

iii = 0;
for i = 1:length(accumulations)
    pos = accumulations(i).positions;
    acum_interpolations.pos_int(iii+(1:length(pos))) = uint32(floor(pos)-1);
    acum_interpolations.pos_frac(iii+(1:length(pos))) = pos-floor(pos);
    %if weights are correctly defined, they are used...
    if(length(accumulations(i).positions) == length(accumulations(i).weights));
        acum_interpolations.weights(iii+(1:length(pos))) = accumulations(i).weights;
    else % by default weights are 1
        acum_interpolations.weights(iii+(1:length(pos))) = ones(length(pos),1);
    end        
        
    acum_interpolations.accum_length(i) = length(pos);
%    acum_interpolations.global_weight(i) = 1/length(pos);
    acum_interpolations.global_weight(i) = accumulations(i).g_weight;
    iii = iii + length(pos);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq_relative alphas] = define_warps_linear_chirps(alpha_max,num_warps,t_orig)
% function [freq_relative alphas] = define_warps_linear_chirps(alpha_max,num_warps,t_orig)
% 
% define warps as relative frequency deviation from original frequency
%
%  in: 
%          alpha_max : maximum value of normalized frequency deviation
%          num_warps : number of warpings
%             t_orig : time vector
%
% out: 
%      freq_relative : relative frequency deviations
%             alphas : chirp rates of the linear chirps
if(1-t_orig(end)*alpha_max <= 0)
    alpha_max = 1/t_orig(end);
    disp('Error: alpha_max and/or nsamps_twarp are too big.');
    disp(['      ' sprintf(' Setting alpha_max equal to %f to continue.',alpha_max)]);
end
alphas = linspace(-alpha_max,alpha_max,num_warps);
freq_relative = zeros(length(t_orig),num_warps);

for i = 1:num_warps
    freq_relative(:,i) = 1+t_orig*alphas(i);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = design_pole(theta,Q,NFFT,db,softness,max_frac_pi)
% function p = design_pole(theta,Q,NFFT,db,softness,max_frac_pi)
%
% Function that designs a pole for the IIR-LTV-Q-FFT transform.
% The Q factor sets the number of cycles of a sinusoid of digital frequency
% theta that must fit into the window. Thus the window is designed so that
% Q cycles fit within its db drop points (normally 3dB). 

if nargin < 4; db = 6; end;
if nargin < 5; softness = 0.01; end;
if nargin < 6; max_frac_pi = 0.6; end;

tau = Q * (2*pi^2) / NFFT / theta;
    
% The hyperbola solutions are calculated
rr = roots([1,-(tau+max_frac_pi*pi),pi*tau*max_frac_pi-softness]);
% the smallest solution is the desired value
tau_h = min(rr);
tau = tau_h;

% The pole is calculated for the saturated Q
w = 1/10^(db/20);
pol = [2*w-(1+cos(tau)) -(4*w*cos(tau)-2*(1+cos(tau))) 2*w-(1+cos(tau))];
r = roots(pol);
% The pole inside the unit circle is the thesired solution
p = r(abs(r)<1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function poles = calculate_poles(thetas,params,NFFT,slope)
% function poles = calculate_poles(thetas,Q,NFFT,slope)
%
% Function that calculates the poles of a IIR-LTV-Q-FFT transform.
% For each digital frequency theta a pole is designed that complies with
% the specified Q factor, for the given number of FFT point.
%
% thetas - digital frequency grid (only half the spectrum - size: NFFT/2 + padding)
%   params.Q - quality factor
%   params.db - db drop
%   params.softness - saturation softness
%   params.max_frac_pi - maximum fraction of window width
%   NFFT - number of points of the FFT


if nargin < 4; slope = 0.9/3; end;

N = length(thetas);
poles = zeros(N,1);

for i=1:N
   poles(i) = design_pole(thetas(i),params.Q+slope*thetas(i),NFFT,params.db,params.softness,params.max_frac_pi);
   %poles(i) = design_pole(thetas(i),Q+slope*thetas(i),NFFT,params);
end

end

