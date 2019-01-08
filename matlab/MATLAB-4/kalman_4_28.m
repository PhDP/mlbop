%-----------------------------------------------------------------
%  MATLAB code for Exercise 4.28 
%  An AR estimation task based on Kalman filter.
%-----------------------------------------------------------------




clear;

rseed	= 1;
rand('state',rseed);

%Create 500 samples of the AR sequence
N = 500;
x = zeros(N,1);
sigma_n = 0.7;
noise = sigma_n*randn(N,1);
a(1) = 0.2; % 0.95   0.8   0.6
a(2) = 0.1;  % 0.9   0.6   0.4
x(1)=4;
x(2)=-2;
for n=3:N
    x(n) = -a(1)*x(n-1) -a(2)*x(n-2) + noise(n);
end;


%Create signal y = x + noise
sigma_v = 0.2;
noise_v = sigma_v*randn(N,1);
y = x + noise_v;


%Kalman Filter
F = [-a(1)  -a(2); 1 0];
H = [1 0];
Q = [sigma_n^2 0; 0 0];
R = sigma_v^2;
%initalization
x_hat = zeros(N,1);
x0 = [0; 0];
%P = ([x(2); x(1)]-x0)*([x(2); x(1)]-x0)';
P = 0.01*[1  0;  0  1];
x_hat(2:-1:1) = x0;
%Iterations
for n=1:N-1
    S = R + H*P*H';
    K = P*H'*inv(S);
    x0 = x0 + K*(y(n) - H*x0);
    P = P - K*H*P;
    x0 = F*x0;
    P = F*P*F' + Q;
    x_hat(n+1:-1:n) = x0;
end;
%end of Kalman Filter


%plot the results
figure(1),plot(x(1:50),'r');hold on; plot(x_hat(1:50),'k');
xlabel('$n$');
ylabel('$x_n$');

figure(2); plot( y(1:50), 'r');
xlabel('$n$');
ylabel('$y_n$');