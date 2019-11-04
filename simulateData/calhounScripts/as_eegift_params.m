%% see the original example at 
%% /export/apps/linux-x86/matlab/toolboxes/GroupICATv2.0e/icatb/toolbox/eegiftv1.0c/Input_eeg_data_subjects_1.m

modalityType = 'eeg';
which_analysis = 2; %% 1- Regular, 2- Group with icasso

if which_analysis==2;
icasso_opts.sel_mode = 'both';  % Options are 'randinit', 'bootstrap' and 'both'
icasso_opts.num_ica_runs = 15; % Number of times ICA will be run
end;

%% Group PCA performance settings. Best setting for each option will be selected based on variable MAX_AVAILABLE_RAM in icatb_defaults.m. 
perfType = 1; % 1- maximize performance, 2 - less memory, 3 - manually spec

dataSelectionMethod = 4; % options 1,2,3,4

if dataSelectionMethod==4;
load params.mat datadir outputDir;

cd(datadir);
files=dir('Subject_*.mat');
for k = 1:length(files);
eval(sprintf('input_data_file_patterns{k} = [''%s'' ''/Subject_%s.mat''];',datadir,num2str(k)));
end;end;
input_data_file_patterns = input_data_file_patterns';


prefix = 'EEGIFT';

%% Group PCA Type. Used for analysis on multiple subjects and sessions.
% Options are 'subject specific' and 'grand mean'. 
%   a. Subject specific - Individual PCA is done on each data-set before group
%   PCA is done.
%   b. Grand Mean - PCA is done on the mean over all data-sets. Each data-set is
%   projected on to the eigen space of the mean before doing group PCA.
%
% NOTE: Grand mean implemented is from FSL Melodic. Make sure that there are
% equal no. of electrodes between data-sets.
%
group_pca_type = 'subject specific'; 

backReconType = 'gica'; % ARBITRARY options 'str' or 'gica'

preproc_type = 4; % ARBITRARY 1-remove mean per tp, 2-remove mean per voxel, 
                  % 3-Intensity normalization, 4-variance normalization,
                  % 5-skip (do it on your own before EEGIFT, NO LONGER
                  % WORKS)

pcaType = 1; % 1- standard, 2- expectation maximization

%% PCA options (Standard)
if pcaType==1;
pca_opts.stack_data = 'yes'; % 'yes'- stacked, lots of memory
                             % 'no'- a pair are loaded at a time
pca_opts.storage = 'full'; % 'full'-full storage of cov matrix, 'packed'-lower triangular

pca_opts.precision = 'double'; %'double' or 'single'

pca_opts.eig_solver = 'selective'; %'selective'
                                   %'all'- all eigs are computed, runs
                                   %%slow, good with convergence issues 
end;

%% ESTIMATION DID NOT WORK
doEstimation = 0;                     %% 1 means estimation of components takes place and that is used
if doEstimation==1; 
estimation_opts.PC1 = 'mean';         %% options are 'mean', 'median', and 'max' for each reduction
estimation_opts.PC2 = 'mean';         %% the length of the cell is equal to # of data reductions
estimation_opts.PC3 = 'mean';
else;end;

numReductionSteps = 2; % max is 3
numOfPC1 = 2; % this is overridden in ss5_eegift.m
numOfPC2 = 2; % this is overridden in ss5_eegift.m
numOfPC3 = 2; % this is overridden in ss5_eegift.m

scaleType = 0; %0- don't scale, 1- scale to data, 2-z scores, 3-normalize using maximum value (and multiple topos using max val)
               %4 - scale topographies using the max value and timecourses using sd of tc

%% algoType = 1; %1-infomax, 2-fastICA, etc. %% this variable is loaded in