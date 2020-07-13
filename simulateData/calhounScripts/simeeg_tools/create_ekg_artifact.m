function [thedata, badChan] = create_ekg_artifact(goodchan, chan_source, chan_sink,sr, a_pwav, lengthK)
%% This function will be used to create EKG artifact in an EEG simulated cap
%% Inputs include:

%% Outputs include the simulated electrode data given the input nodes

%% Create simulated EKG artifact here
% Declare length of simulated EKG here
srm = 1/sr;
x=0:srm:lengthK; % sample frequency
x=0:srm:lengthK;
rate=60; % bpm
li=30/rate;  
%% p wave specifications
a_pwav=a_pwav; % amplitude p wave
d_pwav=0.09; % duration p wave 
t_pwav=0.16; %% inter wave lat
%% q wave specifications
a_qwav=a_pwav*.1;
d_qwav=0.066;
t_qwav=0.166;
%% qrs specifications
a_qrswav=1.6;
d_qrswav=0.11;
%% s wave spec
a_swav=0.25;
d_swav=0.066;
t_swav=0.09;
%% t wav spec
a_twav=0.35;
d_twav=0.142;
t_twav=0.2;
%% u wav spec
a_uwav=a_twav*.01;
d_uwav=0.0476;
t_uwav=0.433;

%% Now create the waves here
pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);
%qwav output
qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);
%qrswav output
qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);
%swav output
swav=s_wav(x,a_swav,d_swav,t_swav,li); 
%twav output
twav=t_wav(x,a_twav,d_twav,t_twav,li);
%uwav output
uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);
%ecg output
ecg=pwav+qrswav+twav+swav+qwav+uwav;


%% Now turn it into EEG data
EEG.data = zeros(64, length(ecg));

for k = 1:length(chan_source)
EEG.data(chan_source(k),:) = ecg./length(chan_source);
end;

for k = 1:length(chan_sink)
EEG.data(chan_sink(k),:) = ecg*-1./length(chan_sink);
end;

EEG = create_eeglab(EEG,goodchan,sr);


%% interpolate
EEG = eeg_interp(EEG,[setdiff(1:64,[chan_source chan_sink])]);

%% Now zero out all negative patterns
% This is a hilarious hockey joke
neg_explore = EEG.data;
neg_explore = neg_explore.';
[valCheck_row valCheck_col] = find(mean(neg_explore)<0);

thedata = EEG.data(EEG.goodchan,:)-repmat(mean(EEG.data(EEG.goodchan,:)),length(EEG.goodchan),1);
badChan = valCheck_col;
end

