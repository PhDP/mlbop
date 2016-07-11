%-----------------------------------------------------------------------------------------------
% Exercise 19.11
% ICA.
% The fast ICA toolbox need to be downloaded from http://research.ics.aalto.fi/ica/fastica/
% and added in the matlab path
% addpath([pwd,filesep,'FastICA']);
%------------------------------------------------------------------------------------------------



function Matlab19_11

 

[x1,fs,nbits] = wavread('voice1.wav');
[x2,fs,nbits] = wavread('voice2.wav');
[x3,fs,nbits] = wavread('music.wav');


X = [x1';x2';x3'];

% Play Original Audio Signals
for i=1:3
   sig = X(i,:)/max(abs(X(i,:) ));
   wavplay(sig,fs); 
end


% Make the mixtures
A=abs(rand(3,3));
Y=A*X;

% Play Audio Mixtures
for i=1:3
   sig = Y(i,:)/max(abs(Y(i,:) ));
   wavplay(sig,fs); 
end

% Perform ICA
[icasig,A,W] = fastica (Y);
G=icasig*X';
[~,Gmax] = max(abs(G));

% Play ICA-based recovered Signals
for i=1:3
    sig = icasig(i,:)/max(abs(icasig(i,:) ));
   wavplay(sig,fs); 
end
    

% Unmix using PCA
[U,S,V] =  svd(Y*Y');
PCAestim = U'*Y;
G = PCAestim * X';

% Play PCA-based recovered Signals
for i=1:3
   sig = PCAestim(i,:)/max(abs(PCAestim(i,:) ));
   wavplay(sig,fs1); 
end
