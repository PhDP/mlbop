%-----------------------------------------------------------------
%  Exercise 12.15
%  Linear prediction with Gaussian prior
%  Use Matlabs errorbar function.
%-----------------------------------------------------------------


clc; clear; close all; format compact; format long eng;

rng('default');

% true signal curve
x = (0:.1e-4:2).'; 
y = .2 * ones(length(x),1) - x + .9 * x.^2 + .7 * x.^3 - .2 * x.^5; 

% plot the true curve
% figure; plot(x1,y1,'b'); 

% training samples (20 or 500)
N = 20; 

% sample interval [a b]
a = 0; b = 2; 

% generate samples
x1 = (a : b/N : b - b/N).';

% noise generation 
sigma_n = .05;
n = sqrt(sigma_n) .* randn(N,1); 

% use the true theta 
theta_true = [.2; -1; .9; .7; -.2]; 

% or a random one
theta_dstrbd = [-0.004; -10.54; 0.465; 0.087; -.093];
l = length(theta_true); 

% compute the measurement matrix
Phi = [ones(N,1) x1 x1.^2 x1.^3 x1.^5];

% generate noisy observations using the linear model
y1 = Phi * theta_true + n;

% set the parameters of Gaussian prior
sigma_theta = 2;
mu_theta_prior = theta_true; % or mu_theta_prior = theta_dstrbd;

% compute the covariance matrix of the Gaussian posterior
Sigma_theta_pos = inv(sigma_theta^-1 * eye(l) + sigma_n^-1 * (Phi.' * Phi));

% compute the posterior mean
mu_theta_pos =  mu_theta_prior + sigma_n^-1 * Sigma_theta_pos * Phi.' * (y1 - Phi * mu_theta_prior); 

% linear prediction
Np = 20; 

% generate prediction samples
x2 = (b-a) * rand(Np,1); 

% compute prediction measurement matrix 
Phip = [ones(Np,1) x2 x2.^2 x2.^3 x2.^5];

% compute the predicted mean and variance
mu_y_pred = Phip * mu_theta_pos;
sigma_y_pred = diag(sigma_n + sigma_n * sigma_theta * Phip * inv(sigma_n * eye(l) + sigma_theta * (Phi.' * Phi) ) * Phip.');

% plot the predictions along the condifence intervals
figure; 
set(gca,'FontSize',12,'FontName','Times');
plot(x,y,'k',x2,mu_y_pred,'kx'); 
hold on; errorbar(x2,mu_y_pred,sigma_y_pred,'r.'); hold off; axis tight; 
xlabel('x'); ylabel('y');
