function topohead(data,varargin);
%% example plot of a single electrode:
%% topohead(zeros(1,30),'sigelectrode',[7 8 18],'goodchan',[],'colormap',[1 1 1],'dotsize',16,'dotcol','k');

%% parse arguments (set the defaults here)
searchstr = {'colormap','maplimits','sigelectrode','style','goodchan','dotsize','dotcol','cap'};
defvals   = [{'jet'}     ,{'absmax'}   ,{[]}            ,{'both'} ,{1:length(data)},{12},{'k'},'Olin'];

%%%%%%%%%% you shouldn't need to change this part
if nargin==1;varargin={'','','','','','','','','',''};else;
for k=(length(varargin)+1):10;varargin{k}='';end;
end;

for k=1:length(searchstr);
thein = find(ismember(varargin(1:2:end),searchstr(k))==1);     
if thein>0; theval{k}=varargin{thein*2};else;theval{k}=defvals{k};end;
end;

if length(data)==64 & strcmp(theval{8},'Olin')==1;
EEG.chanlocs=loadchan_64;
%sel_chan = [1 3 8 12 26 30 46 50 61 63]; %1-Fp1,3-Fp2,8-F3,12-F4,26-C3,30-C4,46-P3,50-P4,61-O1,63-O2
sel_chan = []; %% disabled for my simulation
else;end;
%%%%%%%%%%

if length(theval{3})>0;sel_chan=theval{3};theval{1}=ones(64,3);else;end;

topoplot(data,EEG.chanlocs,'colormap',theval{1},'maplimits',theval{2},'emarker',{'.',theval{7},theval{6},1},'emarker2',{sel_chan,'.','k',24,1},'style',theval{4},'plotchans',theval{5});
set(gcf,'Color','w');