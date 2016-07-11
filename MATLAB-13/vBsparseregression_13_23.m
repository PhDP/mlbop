%-----------------------------------------------------------------
%  Exercise 13.23
%  Variational Bayes Sparse Regression 
%-----------------------------------------------------------------


clc; clear; close all; format long eng; format compact;

rng('default');

% #sampling points
M = 100; 

% SNR value
SNR = 6; 

% positions of the non-zero basis
K = [22; 64]; 

% generate sampling points
x = linspace(-10,10,M).'; 

% basis matrix
Phi = zeros(M);
for i = 1 : M
    Phi(:,i) = exp(- .5 * 10 * (x-x(i)).^2);
end
Phi_cov = Phi.' * Phi;

% original signal
y_0 = Phi(:,K(1)) + Phi(:,K(2)); 

% add noise
s_y = y_0.' * y_0 /M;
s_n = s_y / 10^(SNR/10);
n = sqrt(s_n) * randn(M,1);

% observations
y = y_0 + n; 

% least squares estimation
w_LS = (Phi_cov) \ (Phi.' * y);
y_LS = Phi * w_LS;

% EM algorithm
% initialization
EMiter = 100;
beta = 1e0;
alpha_EM = 1;
Sigma_EM = eye(M);
for i = 1 : EMiter
    Sigma_EM = inv(beta * (Phi_cov) + alpha_EM * eye(M));
    mu_EM = beta * Sigma_EM * Phi.' * y;
    
    alpha_EM = M/(mu_EM.' * mu_EM + trace(Sigma_EM));
    
    beta = M / (norm(y - Phi * mu_EM)^2 + trace(Sigma_EM * Phi_cov));
end
% EM estimates
y_EM = Phi * mu_EM; 


% Variational Bayes algorithm
% initialization
VBiter = 200;
beta = 1e2;
alpha_VB = 3e5 * ones(M,1);
hyp_j = 0;
mu_VB = zeros(M,1);
for i = 1 : VBiter
    Sigma_VB = inv(beta * (Phi_cov) + diag(alpha_VB));   
    mu_VB = beta * Sigma_VB * Phi.' * y;
    
    alpha_VB = (hyp_j + .5)./(hyp_j + .5 * (mu_VB.^2 + diag(Sigma_VB)) );
    
    beta = (hyp_j + .5 * M)/(hyp_j + .5 * (norm(y - Phi * mu_VB)^2 + trace( Sigma_VB * Phi_cov)) );
end
% variational Bayes estimates
y_VB = Phi * mu_VB; 

% plot the results
figure; 
xx = -10:.01:10;
% use spline to smooth curves
y_0p = spline(x,y_0,xx); 
y_LSp = spline(x,y_LS,xx);
y_EMp = spline(x,y_EM,xx);
y_VBp = spline(x,y_VB,xx);
plot(x,y,'xk'); hold on; plot(xx,y_0p,'r',xx,y_LSp,'--r',xx,y_EMp,'--k',xx,y_VBp,'k'); hold off; legend('Observations','Original','ML','EM','Variational');
