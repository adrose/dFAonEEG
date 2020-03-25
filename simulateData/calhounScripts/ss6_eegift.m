function ss6_eegift(H,tg,tag);
%% Identify the EEGIFT source which best matches with simulated source and display

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


map=CMRmap;


%% obtain the frequencies and ground truth (for flipping)
cd(thepath);eval(sprintf('cd(''ss3_plot_sources/ss3_plot_sources_dur%d_olap%d/Sub_1'');',dur,olap));
load ss3_plot_sources_sub01.mat sources_sim_fft thefreq goodchan dur

[chans bins tt nums]=size(sources_sim_fft);

sources1_gt = squeeze(mean(sources_sim_fft,3));


outputDir = [thepath 'OUTPUT_' tag]; 

cd(outputDir);
load EEGIFT_ica.mat W icasig mask_ind A


%% create a sim sources vs ica sources correlation matrix (corrmat)
for G = 1:size(icasig,1)
    themat = reshape(icasig(G,:),bins,chans)';
    
    
    %% compute correlation with the ground truth
    thevec1 = reshape(themat,bins*chans,1);
    for k = 1:(nums-1)
    thevec2 = reshape(sources1_gt(:,:,k),bins*chans,1);
    corrmat(G,k)=corr(thevec1,thevec2); 
    end;
end;

%% determine which ica source best matches which sim source (gtmatch)
for k = 1:(nums-1)
[theval thein]=max(abs(corrmat(:,k)));
gtmatch(k) = thein;               %% keep track of which gt matches
end;

%% reshape the ICA components and determine whether they should be flipped
for G = 1:size(icasig,1)
[theval2 thein2]=max(abs(corrmat(G,:)));    
fliplist(G) = sign(corrmat(G,thein2));

ICA_comps(:,:,G) = reshape(icasig(G,:),bins,chans)';   
ICA_comps(:,:,G) = ICA_comps(:,:,G).*fliplist(G);
end;


% %% plot the correlation with each source for each component 
% FH=figure;set(FH,'visible','off');
% for j = 1:size(corrmat,1)
% for k = 1:size(corrmat,2)
% eval(sprintf('text(%d,%d,''%s'',''fontweight'',''bold'');',k,corrmat(j,k),num2str(j)));hold on;
% %eval(sprintf('text(k,abs(corrmat(j,k)),''%s'',''fontweight'',''bold'');',num2str(j)));hold on;
% end;end;
% ylim([-.05 1]);xlim([.5 4.5]);
% 
% set(FH,'Color','w');
% cd(outputDir);
% eval(sprintf('export_fig Correlation_Sim_and_Comp_dur%d.png',dur));


for G = 1:size(icasig,1);
    
    data = zeros(64,bins);
    data(goodchan,:)=ICA_comps(:,:,G);
    
%     %% plot the spectrum
%     FH=figure;
%     plot(thefreq,mean(data(goodchan,:)),'k','linewidth',2);
% 
%     axis square;
%     set(FH,'Color','w');
%     set(gca,'xtick',[0 4 8 16 32]);
%     cd(outputDir);
%     eval(sprintf('export_fig Source_Spectrum_%d_dur%d.png -m2.5',G,dur));
% 
    %% plot the average topography
    FH=figure;set(FH,'visible','off');
    topohead(mean(data,2),'goodchan',goodchan,'cap','Olin');
    colormap(map);

    cd(outputDir);
    eval(sprintf('export_fig Source_Topography_%d_dur%d.png -m2.5',G,dur));

    %% plot the map
    FH=figure;set(FH,'visible','off');
    
    therange = [min(min(data)) max(max(data))];
    
    chanfreqplot(data',goodchan,therange,thefreq,'olinrest');
    colormap(map);
    
    cd(outputDir);
    eval(sprintf('export_fig Source_Freq_chan_%d_dur%d.png -m2.5',G,dur));
    close all;

end;

save ss6_eegift.mat ICA_comps fliplist corrmat gtmatch

cd(thepath);