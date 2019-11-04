function [h1 newdata chanorder]=chanfreqplot(data,goodchan,range,freq,env,tosmooth);
%% CHANFREQPLOT plots frequency X channels data as a 2D image.
%% Inputs:
%% data: is a freq X channels matrix
%% range: 1X2 upper and lower scale values
%% spikes: solid lines to separate channels

if isempty(range)==1;minx=min(min(data));maxx=max(max(data));
    range=[-max(abs([minx maxx])) (max(abs([minx maxx])))];end;
if nargin<6;tosmooth=0;end;
if nargin<5;env='fmrieeg';end;
if strcmp(env,'fmrieeg');
    thelablist = {'O','P','T','C','F'};
    locs{1} = [9 20 10];                       %% OCCIPITAL
    locs{2} = [16 8 19 7 15];                  %% PARIETAL
    locs{3} = [29 13 14 30];                   %% TEMPORAL
    locs{4} = [28 6 24 18 23 5 27];            %% CENTRAL
    locs{5} = [25 11 21 3 1 17 22 4 2 12 26];  %% FRONTAL
elseif strcmp(env,'olinrest');
    thelablist = {'O','P/C','T','F'};
    locs{1} = intersect([64:-1:53],goodchan);                                 %% OCCIPITAL
    locs{2} = intersect([52:-1:44 40 39 38 37 36 30 29 28 27 26],goodchan);   %% PARIETAL/CENTRAL
    locs{3} = intersect([43 42 41 35 34 33 32 31 25 24],goodchan);            %% TEMPORAL
    locs{4} = intersect([23:-1:1],goodchan);                                  %% FRONTAL
end; 

freqmark = [4 8 13 30 max(freq)];

freqmark = fliplr(freqmark); %% used to put lowest frequencies on bottom
freq = fliplr(freq);         %% used to put lowest frequencies on bottom
lw = 2;

newdata=[];chanorder=[];spike=zeros(1,length(locs));count=0;
for k = 1:length(locs); 
newdata = [newdata data(:,locs{k})]; 
chanorder = [chanorder locs{k}];
count=count+length(locs{k});
spike(k) = count;
end;

newdata = flipud(newdata); %% used to put lowest frequencies on bottom

if tosmooth==0;
    h1=imagesc(newdata,range);hold on;
elseif tosmooth==1;
    newdata=smooth(newdata);
    h1=imagesc(newdata,range);hold on;
else;end;
    
set(gca,'linewidth',2);

%% plot lines to distinguish the channels
spike = [1 spike];
for k = 1:length(spike)-1;
plot([spike(k+1)+.5 spike(k+1)+.5],[0 size(data,1)+.5],'k','linewidth',lw);
thelabloc(k) = sum(spike(k)+spike(k+1))/2;
end;

% %% new code to plot lines to distinguish frequencies
for k = 1:length(freqmark);
[tempval temppos] = min(abs(freqmark(k)-freq)); 
theloc(k)=temppos;
plot([0 size(data,2)+.5],[theloc(k) theloc(k)],'k','linewidth',lw);
end;
set(gca,'Ytick',theloc,'Yticklabel',freqmark);

set(gca,'Xtick',thelabloc,'Xticklabel',thelablist);
ylim([0 size(newdata,1)]+.5);xlim([.5 size(newdata,2)+.5]);
set(gca,'fontsize',12);
