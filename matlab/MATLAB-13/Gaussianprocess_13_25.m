%-----------------------------------------------------------------
%  Exercise 13.25
%  Gaussian process regression
%-----------------------------------------------------------------


clc; clear; close all; format long eng; format compact;

rng('default');

% generate samples from a GP first
N = 20; 
% input samples
x =  randn(N,1); 
% plot samples
% plot(x,ones(1,N),'.');

% compute the gaussian kernel
[x1,x2] = meshgrid(x);

% length scale parameter
h = .5;
K = exp(- .5 * (x1-x2).^2 ./ h^2 ) ; 

% noise variance
noise_var = 0.1;

% generate data 
y = chol(K).' * randn(N,1) + sqrt(noise_var) * randn(N,1);
figure; plot(x, y, '+')

% prediction points
D = 1e2;
xp = linspace(-3, 4, D).'; 

% Compute the mean and the covariance function of the predictive GP
Sigma_N = K + noise_var * eye(N);
xk1 = repmat(x, [1, D]);
xk2 = repmat(xp.', [N, 1]);
k = exp(- .5 * (xk1 - xk2).^2 ./ h^2 ) ;
muf = zeros(D,1);
sigma_y = zeros(D,1);
Ry = Sigma_N \ y;
for d = 1 : D
    muf(d) = k(:,d).' * Ry;
    sigma_y(d) =  1 - k(:,d).' * ( Sigma_N \ k(:,d)); 
end
f = [muf + 2 * sqrt(sigma_y); flipdim(muf-2*sqrt(sigma_y),1)];

% plot the results 
fill([xp; flipdim(xp,1)], f, [7 7 7]/8);
figure; hold on; plot(xp, muf, 'r'); plot(x, y, '+k'), axis tight;

% 2nd version - Compute the mean and the covariance function 
% of the predictive GP using matrix calculations - 
[xxp1,xxp2] = meshgrid(x,xp);
Kxxp = exp(- .5 * (xxp1-xxp2).^2 ./ h^2 );
[xp1,xp2] = meshgrid(xp,xp); 
Kxp = exp(- .5 * (xp1-xp2).^2 ./ h^2 ); 

K_f = [K+noise_var*eye(N) Kxxp.'; Kxxp Kxp];
mu_f = Kxxp * inv(K+noise_var*eye(N)) * y;
Sigma_f = Kxp - Kxxp * inv(K+noise_var*eye(N)) * Kxxp.';
f = [mu_f + 2 * sqrt(diag(Sigma_f)); flipdim(mu_f-2*sqrt(diag(Sigma_f)),1)];

% plot the results 
figure;
fill([xp; flipdim(xp,1)], f, [7 7 7]/8);
hold on; plot(xp, mu_f, 'r'); plot(x, y, '+k'), axis tight;
