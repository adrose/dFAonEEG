function [thedata tempC_sim] = simulated_eeg(themu,thesigma,wav_type,sr,tempL1,coeff_weight,goodchan,chan_source,chan_sink);

%% create a simulated signal
thein = [0; cumsum(tempL1)];
tempC_sim = [];

for k = 1:length(themu)
into{k} = [(thein(k)+1):thein(k+1)];     
tempC_sim = [tempC_sim random('logistic',themu(k),thesigma(k),[1 length(into{k})]).*coeff_weight(k)];
end;

%%%%%%%%%% RECONSTRUCT SIMULATED
%% 'd' is detail, 'a' is approximation, and Vars.level is the space you want.
recon_sim(1,:) = wrcoef('a',tempC_sim,tempL1,wav_type,4); %delta (0-3.75 Hz)
recon_sim(2,:) = wrcoef('d',tempC_sim,tempL1,wav_type,4); %theta (3.75-7.5)
recon_sim(3,:) = wrcoef('d',tempC_sim,tempL1,wav_type,3); %alpha (7.5-15) "this encompasis the alpha range approximately"
recon_sim(4,:) = wrcoef('d',tempC_sim,tempL1,wav_type,2); %beta (15-32

tempsig = wrcoef('d',tempC_sim,tempL1,wav_type,1);
[tempC3 tempL3] = wavedec(tempsig,1,wav_type); 
recon_sim(5,:) = wrcoef('a',tempC3,tempL3,wav_type,1); %gamma (32-48)ish 

recon_sim_sig = sum(recon_sim);


%% create the format that makes EEGLAB happy
EEG.data = zeros(64,length(recon_sim_sig));


for k = 1:length(chan_source)
EEG.data(chan_source(k),:) = recon_sim_sig./length(chan_source);
end;

for k = 1:length(chan_sink)
EEG.data(chan_sink(k),:) = recon_sim_sig*-1./length(chan_sink);
end;

EEG = create_eeglab(EEG,goodchan,sr);


%% interpolate
EEG = eeg_interp(EEG,[setdiff(1:64,[chan_source chan_sink])]);

%% average reference (I think it makes sense to average reference individually)
thedata = EEG.data(EEG.goodchan,:)-repmat(mean(EEG.data(EEG.goodchan,:)),length(EEG.goodchan),1);




