%-----------------------------------------------------------------
%  Exercise 4.26
%  Noise cancelation 
%-----------------------------------------------------------------



clear;
rseed	= 1;
rand('state',rseed);


%create 5000 data samples of a simple sinusoidal signal
N = 5000;
omega0 = 2*10^(-3)*pi;
s = cos(omega0*(1:N) )';
%plot(s);


%Create 5000 data samples of the AR proccesses  v_1 and v2
v1 = zeros(N,1);
v2 = zeros(N,1);
d = zeros(N,1);

sigma = 0.05;
noise = sigma*randn(N,1);
a1= 0.8;
a2= 0.5;

v1(1)=0;
v2(1)=0;
for n=2:N
    v1(n) = a1*v1(n-1) + noise(n);
    v2(n) = a2*v2(n-1) + noise(n);
end;

%Solve the Least Squares problem
A = [ sigma^2/(1-a2^2)   a2*sigma^2/(1-a2^2);   a2*sigma^2/(1-a2^2)  sigma^2/(1-a2^2)];
b = [ sigma^2/(1-a1*a2);  a1*sigma^2/(1-a1*a2)];
w= A\b;

%Create the sequence of the "restored" signal
d(1) = w(1)*v2(1);
for n=2:N
    d(n) = w(1)*v2(n) + w(2)*v2(n-1);
end;

% figure(1);
% subplot(1,2,1);
% plot(s+v1);
% subplot(1,2,2);
% plot(s+v1-d);

%plot the results
figure(1);
plot(s+v1,'r');
title(['Signal $d_n=s_n+v_1(n)$'], 'Interpreter','latex', 'fontsize', 16);


figure(2);
plot(s+v1-d,'r');
title(['Signal $s_n+v_1(n)-\hat d(n)$'], 'Interpreter', 'latex', 'fontsize', 16);