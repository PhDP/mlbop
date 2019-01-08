%-----------------------------------------------------------------
%  Exercise 12.16
%  EM linear regression
%  Use Matlabs errorbar function.
%-----------------------------------------------------------------


clc; clear; close all; format compact; format long eng;

rng('default');

% true signal curve
x = (0:.1e-4:2).'; 
y = .2 * ones(length(x),1) - x + .9 * x.^2 + .7 * x.^3 - .2 * x.^5; 

% plot the true curve
% figure; plot(x1,y1,'b'); 

% training samples
N = 500; 

% linear coefficients
K = 5; 

% sample interval [a b]
a = 0; b = 2; 

% generate samples
x1 = (a : b/N : b - b/N).';

% noise generation
sigma_eta = .05; 
n = sqrt(sigma_eta) .* randn(N,1); 

% use true parameter theta
theta_true = [.2; -1; .9; .7; -.2]; 

% compute the measurement matrix
Phi = [ones(N,1) x1 x1.^2 x1.^3 x1.^5];
Phi_gram = Phi.' * Phi;

% generate noisy observations using the linear model
y1 = Phi * theta_true + n;

% EM algorithm 
% initializate parameters
% experiment on different initializations
EMiter = 20; 
betaj = 1;
sigma_eta_EM = ones(EMiter,1);
alphaj = 1;
Phiy = Phi.' * y1;
for i = 1 : EMiter
    Sigma_theta = inv(betaj * Phi_gram + alphaj * eye(K));
    mu_theta = betaj * Sigma_theta * Phiy;
    
    alphaj = K/(norm(mu_theta)^2 + trace(Sigma_theta));
    
    betaj = N / (norm(y1 - Phi * mu_theta)^2 + trace(Sigma_theta * Phi_gram));
    sigma_eta_EM(i) = 1/betaj; 
end

% perform prediction on new samples
Np = 10; 

% generate prediction samples
x2 = (b-a) * rand(Np,1); 

% compute prediction measurement matrix 
Phip = [ones(Np,1) x2 x2.^2 x2.^3 x2.^5];

% compute the predicted mean and variance
y_pred = Phip * mu_theta;
y_pred_var = diag(sigma_eta_EM(end) + sigma_eta_EM(end) * 1/alphaj * Phip * inv(sigma_eta_EM(end) * eye(K) + 1/alphaj * Phi_gram ) * Phip.');

% plot the predictions along the condifence intervals
figure; 
set(gca,'FontSize',12,'FontName','Times');
plot(x,y,'k',x2,y_pred,'kx'); 
hold on; errorbar(x2,y_pred,y_pred_var,'r.'); hold off; axis tight; 
xlabel('x'); ylabel('y');

figure;
plot(1:EMiter,repmat(sigma_eta, [EMiter 1]), 1: EMiter, sigma_eta_EM); axis([1 EMiter 0.048 0.060]); ylabel('Noise variance'); xlabel('Iterations');