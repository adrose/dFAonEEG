function ss5_eegift(H,tg,tag);
%% Run EEGIFT on the simulated data

disp('WARNING:do not run this script multiple times at once, they share a params.m file'); 


%% unpack the variables given to this function in runmultiple1.m
tag=tag{H};
for j=1:length(tg.lab);
label = tg.lab{j};
eval(sprintf('%s=tg.val{%d}(%d);',char(label),j,H));
end;

%% identify the path to save outputs
temppath=which('ss1_wavelet_signal.m');
pathend = max([max(find(temppath=='/')) max(find(temppath=='\'))]); 
thepath=temppath(1:pathend);

datadir = [thepath 'DATA_' tag];
outputDir = [thepath 'OUTPUT_' tag]; 
dircheck(outputDir);%% create directory

cd(thepath);

save params.mat datadir outputDir;

param_file=icatb_setup_analysis([thepath '/as_eegift_params.m']);   
load(param_file);

sesInfo.userInput.algorithm=algoType; 
sesInfo.userInput.scaleType=scaleType;
sesInfo.userInput.preproc_type = preproc_type;
sesInfo.userInput.numComp = numcomp_2;
sesInfo.userInput.numOfPC1 = numOfPC1;
sesInfo.userInput.numOfPC2 = numOfPC2;

if ismember(algoType,1)==0;
    for k = 1:2:28;
        sesInfo.userInput.ICA_Options{k} = sesInfo.userInput.ICA_Options{27};
        sesInfo.userInput.ICA_Options{k+1} = sesInfo.userInput.ICA_Options{28};
    end;
else;
    labin = [{'on'} {'off'}]; 
    sesInfo.userInput.ICA_Options{18} = extended; %% extended ([N] perform tanh() "extended-ICA")
    sesInfo.userInput.ICA_Options{20} = numcomp_2; %% ncomps
    sesInfo.userInput.ICA_Options{24} = labin{sphering}; 
    sesInfo.userInput.ICA_Options{26} = labin{bias};  
end;    


labin = [{'Regular'} {'Spatial-temporal regression'} {'GICA3'} {'GICA'}];
sesInfo.backReconType = labin{backReconType};  

labin = [{'standard'} {'expectation maximization'} {'svd'}];
sesInfo.pcaType = pcaType;

sesInfo = icatb_runAnalysis(sesInfo, 2);       %% Initialize parameters
sesInfo = icatb_runAnalysis(sesInfo, 3);       %% Data reduction
sesInfo = icatb_runAnalysis(sesInfo, 4);       %% Calculate ICA
sesInfo = icatb_runAnalysis(sesInfo, 5);       %% Back-reconstruct
sesInfo = icatb_runAnalysis(sesInfo, 6);       %% Calibrate Components

%% view and save icasso
cd(outputDir);
icasso_file=dir('*_icasso_results.mat');
if length(icasso_file)>0;
load(icasso_file(1).name);
icassoShow(sR, 'L', size(icasig, 1), 'colorlimit', [.8 .9]);
for k=1:5;FH=figure(k);set(FH,'visible','off');eval(sprintf('export_fig(''ICASSO_%02d.jpg'',FH);',k));end;
else;end;

close all;
cd(thepath);