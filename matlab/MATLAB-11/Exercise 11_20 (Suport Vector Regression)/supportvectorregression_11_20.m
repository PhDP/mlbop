
%-----------------------------------------------------------------
% Exercise 11.20
% Support Vector Regression experiment 
% on an audio sequence corrupted by 
% Gaussian noise and outliers
%
% This file uses an SMO implementation to train the SVM plus other
% subroutines for the various kernel methods.
%
% 
% You need a .wav file to run the experiment...
%-----------------------------------------------------------------


clear all;
rng(0);


%--------------------------------------------------------------------
% Reading wav file. x corresponds to time instances (is., x_i in [0,1])
% fs is the sampling frequency
% Replace the name "BladeRunner.wav" in wavread wuth the name of the file
% you intend to use.
%--------------------------------------------------------------------
N=100;
samples = 1000;
indices = 1:samples/N:samples;

start = 100000;
[sound,fs,nbits] = wavread('wavs\BladeRunner.wav',[start,start+samples]);
y=sound(indices,1);
Ts = 1/fs; % h periodos deigmatolipsias
x = [0:samples]'*Ts; % oi xronikes stigmes tis deigmatolipsias
x = x(indices);
%-------------------------------------------------------


%Add white Gaussian noise
snr = 15; %dB
awgn(y,snr);

%add outliers
O = 0.8*max(abs(y));
percent = 0.1;
M = floor(percent*N);
out_ind = randperm(N,M);
outs = sign(randn(M,1))*O;
y(out_ind) = y(out_ind) + outs;


M=length(y);









%---------------------Code for SVR-------------------------------------
C=1;
epsilon=0.003;
kernel_type='gaus';
kernel_params=0.004;
disp('-----------------------');
disp('SMO starting...');
tic;
[a1, a2, b, KKT] = SMO_regression2(x, y, C, epsilon, kernel_type, kernel_params);
toc;
disp('SMO completed');


%Generate regressor
t = [0:samples]'*Ts;
M2 = length(t);
z=zeros(M2,1);
for k=1:M2
    z(k) = 0;
    for l=1:M
        z(k) = z(k) + (a1(l) - a2(l))*kappa(x(l), t(k), kernel_type, kernel_params);
    end;
    z(k) = z(k) + b;
end;


%find support vectors
SV = [];
for n=1:N
    if (abs(a1(n))>0.001)||(abs(a2(n))>0.001)
        SV = [SV n];
    end;
end;
%---end SVR----------------------------------------------------------------







%plot Output
% For SVR
figure(1)
hold on;
%plot(x,y);
xlabel('time in sec');
ylabel('amplitude');

plot(t,z,'r', 'LineWidth', 1);
% plot(t,z+epsilon,'-.','color','k', 'LineWidth', 1.4);
% plot(t,z-epsilon,'-.','color','k', 'LineWidth', 1.4);

plot(x,y,'.','MarkerEdgeColor',0.3*[1 1 1],'MarkerSize',5);
plot(x(SV), y(SV) , 'o', 'Markersize', 6, 'Color', 'k');

title(['SVR output  ', 'C= ', num2str(C,4), ', ', 'Number of SVs: ', num2str(length(SV),3)])    
        
    
    
    
    
    

