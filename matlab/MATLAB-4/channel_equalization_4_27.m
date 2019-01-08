%-----------------------------------------------------------------
%  MATLAB code for Exercise 4.27
%  Channel equalization
%-----------------------------------------------------------------


clear;
rseed	= 1;
rand('state',rseed);


N = 50 + 2;

u = zeros(N,1);
d = zeros(N,1);

%Create a set of 50 equiprobable +-1 samples
s = sign(2*rand(N,1)-1);
sigma = sqrt(var(s));

%create a noise sequence
sigma_n = 1;
%add the noise to the original sequence
noise = sigma_n*randn(N,1);

%Solve the normal equations
S = [1.25*sigma^2 + sigma_n^2   0.5*sigma^2  0;  0.5*sigma^2  1.25*sigma^2 + sigma_n^2  0.5*sigma^2;  0  0.5*sigma^2  1.25*sigma^2 + sigma_n^2];
p = [sigma^2; 0.5*sigma^2; 0];
w = S\p;

%create the reconstructed signal
u(1) = 0;
u(2) = 0;
d(1) = 0;
d(2) = 0;
for n=3:N
    u(n) = 0.5*s(n) + s(n-1) + noise(n);
    d(n-1) = w(1)*u(n) + w(2)*u(n-1) + w(3)*u(n-2);
end;


% figure(1);
% hold on;
% subplot(2,2,1);
% stem(s,'k');
% subplot(2,2,2);
% stem(u,'k');
% subplot(2,2,3);
% stem(d,'k');
% stem(s,'r');



%plot the results
figure(1);
stem(s(2:51),'k');

figure(2)
stem(u(2:51),'k');


errors = sum(abs(sign(d(2:51))-s(2:51)))/2;

valids_p(2:51) = sign(d(2:51))== s(2:51);
errors_p(2:51) = sign(d(2:51))~= s(2:51);

figure(3);
hold on;

for n=2:51
    if sign(d(n)) == s(n) 
        stem(n-1,sign(d(n)),'r');
    else
        stem(n-1,sign(d(n)),'k');
    end;
end;
% stem(sign(d(2:51)),'k');
% stem(s(2:51),'r');

