
%-----------------------------------------------------------------
% Exercise 11.19
% Kernel Ridge Regression  
% on an audio sequence corrupted by 
% Gaussian noise and outliers
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








%----------Code for unbiased L2 Kernel Ridge Regression (KRR-L2)-----------
C=0.0001;
kernel_type='gaus';
kernel_params=0.004;
sol = kernel_regression_l2_l2_unbiased(x', y, C, kernel_type, kernel_params);
a0 = sol(1:N);

%Generate regressor
t = [0:samples]'*Ts;
M2 = length(t);
z0=zeros(N,1);
for k=1:M2
    z0(k) = 0;
    for l=1:N
        z0(k) = z0(k) + a0(l)*kappa(x(l), t(k), kernel_type, kernel_params);
    end;
end;
%------------------end unbiased KRR-L2-------------------------------------






%-----------Code for biased L2 Kernel Ridge Regression (KRR-L2)------------
C=0.0001;
kernel_type='gaus';
kernel_params=0.004;
sol = kernel_regression_l2_l2_biased(x', y, [C], kernel_type, kernel_params);
a1 = sol(1:N);
b1 = sol(N+1);

%Generate regressor
t = [0:samples]'*Ts;
M2 = length(t);
z1=zeros(N,1);
for k=1:M2
    z1(k) = 0;
    for l=1:N
        z1(k) = z1(k) + a1(l)*kappa(x(l), t(k), kernel_type, kernel_params);
    end;
    z1(k) = z1(k) + b1;
end;
%------------------end biased KRR-L2---------------------------------------




% For unbiased KRR-L2
figure(1)
hold on;
%plot(x,y);
xlabel('time in sec');
ylabel('amplitude');

plot(t,z0,'r', 'LineWidth', 1);

plot(x,y,'.','MarkerEdgeColor',0.3*[1 1 1],'MarkerSize',5);

title(['unbiased KRR-L2  ', 'C= ', num2str(C,4)])

    

% For biased KRR-L2
figure(2)
hold on;
%plot(x,y);
xlabel('time in sec');
ylabel('amplitude');

plot(t,z1,'r', 'LineWidth', 1);

plot(x,y,'.','MarkerEdgeColor',0.3*[1 1 1],'MarkerSize',5);

title(['biased KRR-L2  ', 'C= ', num2str(C,4)])
    
    
    
    
    
    
    

