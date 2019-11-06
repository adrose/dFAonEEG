function ss1_wavelet_signal
%% SS1_WAVELET_SIGNAL.m is part of the simeeg toolbox for generating 
%% realistic simulated eeg and conducting spatiospectral group ICA: 
%% http://mialab.mrn.org/software/simeeg

%% SS1_WAVELET_SINGAL.m loads sample eeg, conducts wavelet analysis, and fits 
%% a logistic distribution to the wavelet coefficients within selected frequency
%% bands. 

%% created by David Bridwell 2016_01_14 and used in the following
%% manuscript:
%% Bridwell, D.A., Rachakonda, S., Silva, R.F., Pearlson, G.D., Calhoun, V.D.
%% (2016) Spatiospectral decomposition of multi-subject EEG: evaluating blind 
%% source separation algorithms on real and realistic simulated data.
%% Brain Topography.

%% Use matlab 7.12.0.635 (R2011a) or later

%% Notes:
%% tempC is coefficients, tempC(1:tempL(1)) is A5. The A's are residuals and 
%% the D's are octaves.
%% With 0-128 Hz, you capture 0-64 Hz (Nyquist) and with wavelet you have
%% approximation (A1) and detail (D1) and D1 has 32-64 and A1 has 0-31.
%% (You can't create the entire tree with this code, you can only do it
%% on the approximations). So with 5 levels here you follow the tree down 5
%% from A1 (0-31) to A2 (0-15) all the way to A5 (0-~2). There is a
%% function in matlab which tells you how far you can go. 

%% How to change the mother wavelet? Biorwavf is biorthogonal (type in
%% waveinfo in matlab for general info on mother wavelets, or waveinfo('bior')
%% for info on biorthogonal. We stumbled across 'rbio' as an additional
%% wavelet.

%% Note that dwtmode is a parameter, the default is given by dwtmode in
%% matlab.

%% To use different mother wavelets you change Vars.W, and the [RF,DF] = biorwavf(Vars.W);
%% and [LO_D,HI_D,LO_R,HI_R] =  biorfilt(DF,RF); command    

%% WAVEINFO TELLS YOU HOW TO CALL THE MOTHER WAVELETS AND THE APPROPRIATE
%% FUNCTION. WRCOEF MAY BE A MORE EFFICIENT WAY TO MAKE THAT FUNCTION CALL.

%% You can view the wavelet with this plot(wavefun('bior3.9'),'r');

%% With bior3.9, the 3.9 is changing the degree and the order of the
%% scaling functions which for biorthogonal are spline functions of
%% different degree, waveinfo('bior') tells you the number of entries.

%% A spline takes linear pieces and you can use them to interpolate data.
%% They have nodes and at the nodes you connect linear pieces in such a way
%% that the derivative in connecting the two is the same from both sides
%% and you control how many derivatives are the same. So you control the
%% value at the node for both lines is the same and the derivative is the
%% same, you basically stitch together functions to make them continuous. 


%% identify the path to save outputs
temppath=which('ss1_wavelet_signal.m');
pathend = max([max(find(temppath=='/')) max(find(temppath=='\'))]); 
thepath=temppath(1:pathend);

rand('state',100); 


%%%%%%%%%% Load sample EEG signal. 
%% sample_eeg is channel Fp1 from a single subject with 250 Hz 
%% sampling rate (note that 250 Hz sampling rate is required for the
%% frequency bands listed below to be correct). It was collected during
%% 4 min ~46 seconds of rest. 
load sample_eeg.mat sample_eeg %% 


Vars.W = 'bior3.9'; % bior - biorthogonal spline 
                    % dbN - Daubechies (not biorthogonal), see waveinfo('db'), you call with db1, db2, etc.


Vars.level = 5;      

[tempC1 tempL1] = wavedec(sample_eeg,Vars.level,Vars.W);


%% determine the index to A5, D5, D4, D3, D2, D1 (with Vars.level5)
thein = [0; cumsum(tempL1)];
for k=1:Vars.level+1;
    into{k} = [(thein(k)+1):thein(k+1)];
end;


tempsig = wrcoef('d',tempC1,tempL1,Vars.W,1);
[tempC2 tempL2] = wavedec(tempsig,1,Vars.W); 
recon(5,:) = wrcoef('a',tempC2,tempL2,Vars.W,1); %gamma 32-48ish


for k = 1:(Vars.level+1)
    into{k} = [(thein(k)+1):thein(k+1)];
    coeffs{k} = tempC1(into{k});
    
    %% estimate the parameters of a logistic distribution (superior to normal)
    temppd = fitdist(coeffs{k},'logistic');
    themu(k) = temppd.mu;
    thesigma(k) = temppd.sigma;        
end;


%%%%%%%%%% RECONSTRUCT ORIGINAL
%% 'd' is detail, 'a' is approximation, and Vars.level is the space you want.
recon(1,:) = wrcoef('a',tempC1,tempL1,Vars.W,5); %delta (0-4 Hz) (used to be 4)
recon(2,:) = wrcoef('d',tempC1,tempL1,Vars.W,5); %theta (4-8)    (used to be 4)
recon(3,:) = wrcoef('d',tempC1,tempL1,Vars.W,4); %alpha (8-16) (used to be 3) "this encompasis the alpha range approximately"
recon(4,:) = wrcoef('d',tempC1,tempL1,Vars.W,3); %beta (16-32) (used to be 3)

tempsig = wrcoef('d',tempC1,tempL1,Vars.W,1);
[tempC2 tempL2] = wavedec(tempsig,1,Vars.W); 
recon(5,:) = wrcoef('a',tempC2,tempL2,Vars.W,1); %gamma (32-48)ish

recon_sig = sum(recon(1:4,:));  %% only reconstruct the first 4  

figure;
plot(recon_sig(1:250*2));hold on;
plot(sample_eeg(1:250*2),'r');
legend('wavelet reconstructed','original');


cd(thepath);
save ss1_wavelet_signal.mat



