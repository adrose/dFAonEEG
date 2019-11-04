function EEG = create_eeglab(EEG,goodchan,sr);
%% create the format that makes EEGLAB happy

EEG.datalen = size(EEG.data,2);
EEG.nbchan = size(EEG.data,1);
EEG.goodchan = goodchan;
EEG.rate = sr;
EEG.srate = sr;
chanlocs = loadchan_64;
EEG.chanlocs=chanlocs;
EEG.etc.amica = [];
EEG.icaact = [];
EEG.specdata = [];
EEG.icachansind = [];
EEG.icasphere = [];
EEG.icawinv = [];
EEG.specicaact = [];
EEG.trials=1;
EEG.epoch = 1;
EEG.event = [];
EEG.setname = [];
EEG.icaweights = [];
EEG.xmax = 0; %this value is adjusted in eeglab
EEG.xmin = 0;