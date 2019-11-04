%% Run the functions in the following order:

%% Step 1: Add a path to simeeg_scripts, 
%% e.g. addpath(genpath('E:\simeeg_scripts'));

%% Step 2: Run ss1_wavelet_signal in matlab
%% SS1_WAVELET_SINGAL.m loads sample eeg, conducts wavelet analysis, and fits 
%% a logistic distribution to the wavelet coefficients within selected frequency
%% bands. 

%% Step 3: Run ss2_create_sim in matlab
%% SS2_CREATE_SIM.m creates a simulated signals based on the distribution of 
%% wavelet coefficients derived from real data in ss1_wavelet_signal.m
%% The user can specify the frequency and spatial characteristics of the sources.

%% Step 4: RUN runmultiple1 in matlab to run EEGIFT on the simulated data created above.
%% RUNMULTIPLE1 uses:
%% ss4_mix_sources.m to mix the simulated sources 
%% ss5_eegift.m to run eegift
%% ss6_eegift.m to identify sources which match simulated sources and plot

%% Note that runmultiple1.m is designed to run many eegift runs over many different
%% parameter options, the different options go into different folder names
%% which have the outputs for that particular run. 