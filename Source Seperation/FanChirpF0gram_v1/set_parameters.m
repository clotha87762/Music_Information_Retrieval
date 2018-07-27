% script to set FChT analysis parameters
% see below for a description of each parameter

if(isunix); dir_slash = '/'; else dir_slash = '\'; end    

%% =============  WARPING PARAMETERS  =============
% maximum frequency of interest
fmax = 10000; % Hz
% number of samples of the warped signal frame
warp_params.nsamps_twarp = 2048;
% number of fft points (controls zero-padding)
nfft = warp_params.nsamps_twarp; %nfft = 3*wp.nsamps_twarp/2; 
% warping type
warp_params.warping_type = 'linear';
% maximum value of normalized frequency deviation (alpha)
warp_params.alpha_max = 6;
% number of warpings
warp_params.num_warps = 21;
% oversampling factor
warp_params.fact_over_samp = 2;
% hop in samples
hop = warp_params.fact_over_samp * 256;

%% ============= CONSTANT Q PARAMETERS =============
% whether to compute a constant-Q FChT
cqt_params.use_CQT = 0;
% quality factor Q for the constant-Q FChT
cqt_params.Q = 70;
cqt_params.db = 6;
cqt_params.softness = 0.02;
cqt_params.max_frac_pi = 0.7;

%% =============   F0-GRAM PARAMETERS  =============
% minimun fundamental frequency
f0_params.f0min = 80; % Hz
% number of octaves
f0_params.num_octs = 4;
% number of f0s per octave
f0_params.num_f0s_per_oct = 192;
% number of f0 hypotesis to extract
f0_params.num_f0_hyps = 10;
% type of f0gram returned: 'none', 'single', 'compound'
f0_params.f0gram_type = 'compound';
% whether to save the f0gram
f0_params.save_f0gram = 0;
% whether to use a f0 preference guassian function
f0_params.prefer = 1;
% mean of f0 preference guassian function (MIDI number for C4)
f0_params.prefer_mean = 60;  
% stdev of f0 preference guassian function (stdev in MIDI numbers)
f0_params.prefer_stdev = 18;  

%% ======== GATHERED LOG SPECTRUM PARAMETERS =======
% high-pass logS
glogs_params.HP_logS = 1;
% whether to attenuate subharmonics
glogs_params.att_subharms = 1;
% whether to apply the GlogS correction model 
glogs_params.apply_GlogS_model = 1; % may be changed later
% whether to learn the GlogS correction model from a labeled database
glogs_params.compute_GlogS_model = 0; % may be changed later
% model parameter variables (default values)
glogs_params.median_poly_coefs = [-0.000000058551680 -0.000006945207775 0.002357223226588];
glogs_params.sigma_poly_coefs  = [ 0.000000092782308  0.000057283574898 0.022199903714288];

%% ========  DATABASE AND LABELS PARAMETERS  =======
% whether to use the chirp rate given by F0 labels
% 0 = no, 1 = yes, 2 = calculate it but don't use it
labels_params.chirp_rate_from_labels = 0;
% compute only frames where there is a melody label
labels_params.only_melody_frames = 0;
% labeled database (different databases have different kind of labels)
labels_params.database = 'mirex'; % 'mirex' or 'rwc'
% path to the labeled database 
% set path to audio databases (use set_path_audio_databases script or set it manually here)
% set_path_audio_databases
path_rwc = './';
path_mirex = './';
if strcmp(labels_params.database, 'rwc'); 
    labels_params.database_path = path_rwc; 
else
    labels_params.database_path = path_mirex; 
end
% audio filename (needed to derive labels filename)
labels_params.audio_file = '';
% sampling frequency of files in the database (used in the design)
labels_params.fs_database = 44100;
% path to the .mat files for pitch and f0gram 
labels_params.path_mat = ['data' dir_slash labels_params.database dir_slash];
