function ss4_mix_sources(H,tg,tag);
%% SS4_MIX_SOURCES.m is part of the simeeg toolbox for generating 
%% realistic simulated eeg and conducting spatiospectral group ICA: 
%% http://mialab.mrn.org/software/simeeg

%% It is called within runmultiple1.m

%% SS4_MIX SOURCES.m determines the mixing matrix four sources and mixes them.
%% It generates a dataset that eegift likes, data is [freq X chan] X 1 X epochs
%% for example, data(1:129,1,1) is the spectrum for the first channel for
%% the first epoch. 

%% created by David Bridwell 2016_01_14 and used in the following
%% manuscript:
%% Bridwell, D.A., Rachakonda, S., Silva, R.F., Pearlson, G.D., Calhoun, V.D.
%% (2016) Spatiospectral decomposition of multi-subject EEG: evaluating blind 
%% source separation algorithms on real and realistic simulated data.
%% Brain Topography.

%% This script uses export_fig: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

%% Use matlab 7.12.0.635 (R2011a) or later

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



cd(thepath);cd('ss2_create_sim');
subs=length(dir('ss2_create_sim_sub*.mat')); %% total subjects


for G = 1:subs

 
cd(thepath);
eval(sprintf('cd ''ss3_plot_sources/ss3_plot_sources_dur%d_olap%d/Sub_%d''',dur,olap,G));
eval(sprintf('load ss3_plot_sources_sub%02d.mat sources_sim_fft source_map;',G));

[numchan numfreq epochs sources] = size(sources_sim_fft);

    for J = 1:epochs
            
    %% fit simulated map to source map (determine mixing matrix)
    for H = 1:(sources-1)
    thex = reshape(source_map(:,:,H),numchan*numfreq,1);  
    they = reshape(sources_sim_fft(:,:,J,H),numchan*numfreq,1);
    
    theP = polyfit(thex,they,1);
    themix(J,H)=theP(1); %% the mixing matrix is the multiplier
    end;
     
    %% reconstruct
    recon = 0;
    for H =  1:(sources-1)
    thex = reshape(source_map(:,:,H)',numchan*numfreq,1);  
    recon = recon+themix(J,H)*thex;
    end;

    sources_recon_fft(:,:,J) = reshape(recon,numfreq,numchan)'; %% save as [chanXfreq] X epoch
    data(:,1,J)=recon;                                         %% save for gift
    
    end;
    
%% add noise
themedian = median(reshape(sources_recon_fft,numchan*numfreq*epochs,1));
addnoise1 = randn(size(sources_recon_fft)).*(themedian*thenoise/100);

thevar_sig = var(reshape(sources_recon_fft,numchan*numfreq*epochs,1));
thevar_noise = var(reshape(addnoise1,numchan*numfreq*epochs,1));
thesnr = (thevar_sig+thevar_noise)./thevar_noise;

for k = 1:size(addnoise1,3);
addnoise2(:,1,k) = reshape(addnoise1(:,:,k)',numchan*numfreq,1);
end;

sources_recon_fft = sources_recon_fft+addnoise1;


%% IN CASE YOU SAVE OUT THE OTHER VARIABLES
%cd(thepath);
%eval(sprintf('dircheck(''ss4_mix_sources_dur%d_olap%d_noise%d'');',dur,olap,thenoise));
%eval(sprintf('save ss4_mix_sources_sub%02d.mat themix source_map sources_recon_fft addnoise1 thesnr',G));


%% save data in the format eegift likes
data = data+addnoise2;

cd(thepath);
eval(sprintf('dircheck(''DATA_%s'');',tag));
eval(sprintf('save Subject_%d.mat data',G));

end;

