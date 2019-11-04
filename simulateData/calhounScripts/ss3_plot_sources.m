function ss3_plot_sources(dur,olap)
%% SS3_PLOT_SOURCES.m is part of the simeeg toolbox for generating 
%% realistic simulated eeg and conducting spatiospectral group ICA: 
%% http://mialab.mrn.org/software/simeeg

%% SS3_PLOT_SOURCES.m plots the sources, saving the figures in ss3_plot_sources 
%% (you can find the outputs already in that folder). This script determines
%% ground truth for the mixing in the next script SS4_MIX_SOURCES.m 

%% created by David Bridwell 2016_01_14 and used in the following
%% manuscript:
%% Bridwell, D.A., Rachakonda, S., Silva, R.F., Pearlson, G.D., Calhoun, V.D.
%% (2016) Spatiospectral decomposition of multi-subject EEG: evaluating blind 
%% source separation algorithms on real and realistic simulated data.
%% Brain Topography.

%% This script uses export_fig: http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig

%% Use matlab 7.12.0.635 (R2011a) or later
if nargin<2;
    dur = 750; %% samples per interval (3 seconds for 250 Hz sample rate)
    olap = 75; %% % overlap between sucessive samples
else;end;

%% identify the path to save outputs
temppath=which('ss1_wavelet_signal.m');
pathend = max([max(find(temppath=='/')) max(find(temppath=='\'))]); 
thepath=temppath(1:pathend);


cd(thepath);
addpath(genpath([thepath 'export_fig']));
addpath(genpath([thepath 'simeeg_tools']));
disp('WARNING: SIMEEG ADDED export_fig TO THE PATH--POSSIBLY REPLACING YOUR CURRENT VERSION');


map = CMRmap;

cd ss2_create_sim;
subs=length(dir('ss2_create_sim_sub*.mat')); %% total subjects


for F = 1:subs;


cd(thepath);cd('ss2_create_sim');
eval(sprintf('load ss2_create_sim_sub%02d.mat sources_sim sr goodchan;',F));

mixture_sim = sum(sources_sim,3);
sources_sim(:,:,size(sources_sim,3)+1)=mixture_sim;


onsets = round(1:(dur-olap/100*dur):(size(sources_sim,2)-dur)); %% considers overlap

for G = 1:size(sources_sim,3);
for H = 1:(length(onsets))  

into = onsets(H):(onsets(H)+dur-1);

dataseg = sources_sim(:,into,G);


[ticafft icafreq] = getfft(dataseg,sr,35,2);

sources_sim_fft(:,:,H,G)=log(abs(ticafft));    


end;
end;


thefreq = icafreq;
data = zeros(64,length(thefreq));


for G = 1:size(sources_sim,3);
data(goodchan,:) = mean(sources_sim_fft(:,:,:,G),3);


%% plot the spectrum
FH=figure;
plot(thefreq,mean(data(goodchan,:)),'k','linewidth',2);

ylabel('log amplitude spectrum');
xlabel('Frequency (Hz)');
set(gca,'fontsize',16);

axis square;
set(FH,'Color','w');
set(gca,'xtick',[0 4 8 16 32]);

cd(thepath);cd('ss3_plot_sources');eval(sprintf('dircheck(''ss3_plot_sources_dur%d_olap%d'');',dur,olap));
eval(sprintf('dircheck(''Sub_%s'');',num2str(F)));
eval(sprintf('export_fig Spectrum_%d_dur%d.png',G,dur));


%% plot the average topography
FH=figure;set(FH,'visible','off');
thedat = mean(data,2);

therange2 = [min(thedat(goodchan)) max(thedat(goodchan))];
[vals5 vals6] = simplevals(therange2(1));
[vals7 vals8] = simplevals(therange2(2));

topohead(thedat,'goodchan',goodchan,'cap','Olin','maplimits',therange2);
colormap(map);

set(FH,'Color','w');

cd(thepath);cd('ss3_plot_sources');eval(sprintf('dircheck(''ss3_plot_sources_dur%d_olap%d'');',dur,olap));
eval(sprintf('dircheck(''Sub_%s'');',num2str(F)));
eval(sprintf('export_fig Topography_%d_dur%d_thr1_%d_%d_thr2_%d_%d.png -m2.5',G,dur,vals5,vals6,vals7,vals8));


%% plot the map
FH=figure;set(FH,'visible','off');

therange = [min(min(data)) max(max(data))];
[vals1 vals2] = simplevals(therange(1));
[vals3 vals4] = simplevals(therange(2));

chanfreqplot(data',goodchan,therange,thefreq,'olinrest');
colormap(map);
set(FH,'Color','w');

eval(sprintf('export_fig Freq_chan_%d_dur%d_thr1_%d_%d_thr2_%d_%d.png -m2.5',G,dur,vals1,vals2,vals3,vals4));
close all;


%% Plot a histogram
FH=figure;set(FH,'visible','off');

[siz1 siz2] = size(data(goodchan,:));

    [heights,locations] = hist(reshape(data(goodchan,:),siz1*siz2,1),200);
    width = locations(2)-locations(1);
    heights = heights/(siz1*siz2*width);
    bar(locations,heights,'hist');hold on;
    set(get(gca,'child'),'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
   
    [thef,thexi] = ksdensity(reshape(data(goodchan,:),siz1*siz2,1));
    
    plot(thexi,thef,'k','linewidth',2);
    
    axis square;
    ylim([0 .55]);
    xlim([locations(1) locations(end)]);
    
    thek(G)=kurtosis(reshape(data(goodchan,:),siz1*siz2,1));
      
    set(FH,'Color','w');
    eval(sprintf('export_fig Source_hist_%d_dur%d.png -m2.5',G,dur));


source_map(:,:,G) = data(goodchan,:); %% save the source map, this is the ground truth
end;


% plot the separate spectrums together
FH=figure;
for G = 1:size(sources_sim,3)-1;
plot(thefreq,mean(mean(sources_sim_fft(:,:,:,G),3)),'k','linewidth',2);hold on;

ylabel('log amplitude spectrum');
xlabel('Frequency (Hz)');
set(gca,'fontsize',16);

axis square;
set(FH,'Color','w');
set(gca,'xtick',[0 4 8 16 32]);

eval(sprintf('export_fig Separate_Spectrum_dur%d.png',dur));
end;


% plot the sum of the spectrums (to compare with the spectrum of the sum etc)
FH=figure;
plot(thefreq,mean(mean(sources_sim_fft(:,:,:,end),3)),'k','linewidth',2);hold on;

ylabel('log amplitude spectrum');
xlabel('Frequency (Hz)');
set(gca,'fontsize',16);

axis square;
set(FH,'Color','w');
set(gca,'xtick',[0 4 8 16 32]);

eval(sprintf('export_fig Sum_Spectrum_dur%d.png',dur));
close all;


eval(sprintf('save ss3_plot_sources_sub%02d.mat sources_sim_fft source_map thefreq goodchan dur olap thek',F));


end;

cd(thepath);


