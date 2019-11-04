function ss2_create_sim
%% SS2_CREATE_SIM.m is part of the simeeg toolbox for generating 
%% realistic simulated eeg and conducting spatiospectral group ICA: 
%% http://mialab.mrn.org/software/simeeg

%% SS2_CREATE_SIM.m creates a simulated signals based on the distribution of 
%% wavelet coefficients derived from real data in ss1_wavelet_signal.m
%% The user can specify the frequency and spatial characteristics of the sources 
%% below.

%% created by David Bridwell 2016_01_14 and used in the following
%% manuscript:
%% Bridwell, D.A., Rachakonda, S., Silva, R.F., Pearlson, G.D., Calhoun, V.D.
%% (2016) Spatiospectral decomposition of multi-subject EEG: evaluating blind 
%% source separation algorithms on real and realistic simulated data.
%% Brain Topography.

%% This script uses EEGLAB: http://sccn.ucsd.edu/eeglab/ for interpolation
%% Delorme, A., Makeig, S. (2004) EEGLAB: an open source toolbox for analysis 
%% of single-trial EEG dynamics. Journal of Neuroscience Methods 134:9-21 

%% Use matlab 7.12.0.635 (R2011a) or later

%% identify the path to save outputs
temppath=which('ss1_wavelet_signal.m');
pathend = max([max(find(temppath=='/')) max(find(temppath=='\'))]); 
thepath=temppath(1:pathend);


%% use the version of eeglab within these scripts (i.e. guaranteed to work)
addpath(genpath([thepath 'simeeg_tools']));
addpath(genpath([thepath 'eeglab13_4_4b']))
disp('WARNING: SIMEEG ADDED eeglab13_4_4b TO THE PATH--POSSIBLY REPLACING YOUR CURRENT VERSION');


cd(thepath);
load ss1_wavelet_signal.mat themu thesigma Vars into tempL1;


sr = 250;                               %% sample rate
goodchan = setdiff(1:64,[60 64 33 43]); %% remove I1, I2, M1, M2 to simplify topographic representations

subs = 6;                               %% number of simulated subjects
                            

for H = 1:subs
       
    %% four simulated sources using only the first 4 coefficients, for max frequency of 32 Hz)
    coeff_weight{1} = [1 0 0 0 0 0]; %% scale coefficients
    chan_source{1} = 10;             %% Fz
    chan_sink{1} = 28;               %% Cz

    coeff_weight{2} = [0 1 0 0 0 0]; %% scale coefficients
    chan_source{2} = [7 13];         %% F5 F6
    chan_sink{2} = [54 58];          %% P05 P06

    coeff_weight{3} = [0 0 1 0 0 0]; %% scale coefficients
    chan_source{3} = [58 7];         %% P06 F5
    chan_sink{3} =   28;             %% Cz

    coeff_weight{4} = [0 0 0 1 0 0]; %% scale coefficients
    chan_source{4} = [24];            %% T7
    chan_sink{4} = [32];              %% T8

for G = 1:length(coeff_weight)
[tempsig tempcoeff] = simulated_eeg(themu,thesigma,Vars.W,sr,tempL1,coeff_weight{G},goodchan,chan_source{G},chan_sink{G});
sources_sim(:,:,G) = tempsig;
sources_coeff(:,G) = tempcoeff;
end;

cd(thepath);
cd ss2_create_sim
eval(sprintf('save ss2_create_sim_sub%02d.mat',H));
cd ../

end;







