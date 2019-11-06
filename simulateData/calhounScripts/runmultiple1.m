function runmultiple1;
%% Run on the simulated data.


thisdir=pwd;

eegiftdir = which('eegift');
if length(eegiftdir)==0;
    error('you need to download and add a path for eegift');
else;end;


for numc = [5]; %[5 10 15 20]; %used in the paper
    
for H = 2;%1:2; %1-Infomax, 2-all other algorithms (Infomax is separate since there are additional options)
clear tg

tg.lab{1} = {'dur'};  % duration of window
tg.val{1} = [750];

tg.lab{2} = {'olap'};  % overlap
tg.val{2} = [75];

tg.lab{3} = {'thenoise'};  % variance of the noise is noise % of the median of the 3D data
tg.val{3} = [1]; %[1 5 10 15 20]; %used in the paper

tg.lab{4} = {'numcomp_2'};  % number of components  
tg.val{4} = [numc]; 

tg.lab{5} = {'numOfPC1'};
tg.val{5} = [20]; 

tg.lab{6} = {'numOfPC2'}; %% values with 12 seemed to be unstable (in simulation with 5 sources)
tg.val{6} = [numc]; %used 12 in previous studies (didn't like more than 100)

tg.lab{7} = {'scaleType'}; %% 0 SEEMS TO BE THE BEST (2 IS THE SAME)
tg.val{7} = [0]; %%0-dont' scale, %1-scale to data, 2-z-scores (eegift default) 3-normalize using maximum value (and multiple topos using max val)
               %4-scale topographies using the max value and timecourses using sd of tc

tg.lab{8} = {'preproc_type'}; 
tg.val{8} = [1]; % 1-remove mean per tp, 2-remove mean per voxel, 3-Intensity normalization, 4-variance normalization

tg.lab{9} = {'algoType'}; 
if H==1;
tg.val{9} = [1];  
elseif H==2; 
tg.val{9} = 10; %%  2-'Fast ICA', 3-'Erica', 4-'Simbec', 5-'Evd', 6-'Jade Opac', 7-'Amuse', 8-'SDD ICA', 9-'Radical ICA', 10-'Combi', 11-'ICA-EBM', 12-'FBSS', 13-'IVA-GL'  (see icatb_icaAlgorithm.m and icatb_icaOptions.m)
elseif H==2;         %% (WASOBI is manually done by running with 10, adjusting the call in icatb_icaAlgorithm, and then changing the filenames so it doesnt' write over COMBI)
else;end;

tg.lab{10} = {'backReconType'}; % dont use regular, it makes the individual sum to group, Srinivas says use str or gica
tg.val{10} = [2]; %% 'Regular', 'Spatial-temporal regression', 'GICA3', 'GICA'

tg.lab{11} = {'pcaType'}; %% keep as standard 
tg.val{11} = [1]; %%'standard', 'expectation','svd'

if length(tg.val{9})==1 & tg.val{9}(1)==1;
tg.lab{12} = {'extended'};  
tg.val{12} = [0 1 2]; %[0 1 2];

tg.lab{13} = {'sphering'};  %% keep on
tg.val{13} = [1]; %1-on 2-off

tg.lab{14} = {'bias'};  %% keep on
tg.val{14} = [1]; %1-on 2-off

else;end;


[tg tag] = expandtg(tg);  

for H = 1:length(tag);

%ss4_mix_sources(H,tg,tag);    
%ss5_eegift(H,tg,tag); %% run eegift
ss6_eegift(H,tg,tag); %% visualize components (group) and match to sources

end;
end;
end;




cd(thisdir);
