%-----------------------------------------------------------------
%  Exercise 12.17
%  Gaussian mixture EM algorithm
%  Use error_ellipse function.
%-----------------------------------------------------------------


clc; clear; close all; format long eng; format compact;

rng('default');

% number of samples
N = 300; 

% dimension
l = 2;   

% # of Gaussians
K = 3;  

% # EM iterations 
NofIter = 100; 

% generate data points out of three Gaussians 
% 1st gaussian
N1 = 100; % samples
mu1 = [10 3];
Sigma1 = 1 * eye(l);
R1 = chol(Sigma1);
x1 = repmat(mu1,N1,1) + randn(N1,l) * R1;

% 2nd gaussian
N2 = 100;
mu2 = [1 1];
Sigma2 = 1.5 * eye(l);
R2 = chol(Sigma2);
x2 = repmat(mu2,N2,1) + randn(N2,l) * R2;

% 3rd gaussian
N3 = 100;
mu3 = [5 4];
Sigma3 = 2 * eye(l);
R3 = chol(Sigma3);
x3 = repmat(mu3,N3,1) + randn(N3,l) * R3;

% data point matrix
X = [x1 ; x2; x3].';


% Gaussian mixtures EM algorithm

% initialization
Pk = 1/K* ones(K,1); % initial probabilities

% starting mean, use a random mean 
mu = rand(l,K); 
% or  the true 
% mu = ([3 2 4; 5 .4 3] .* ones(l,K));
% or use a far-off mean 
% mu = [10 11 13; 13 12 11] .* ones(l,K); 
mu = repmat(mu,[ 1 1 NofIter]);

% use an identity covariance matrix as a starting point
Sigma = eye(l); 
Sigmak = repmat(Sigma,[1 1 K]);
Sigmak_inv = zeros(size(Sigmak));
Sigmak_det = zeros(K,1);
gammakn = zeros(K,N);
lg_lklhd = zeros(NofIter,1); % initialize log likelihood variable

% plot data and initial estimates
figure;
plot(x1(:,1),x1(:,2),'.k',x2(:,1),x2(:,2),'+k',x3(:,1),x3(:,2),'*k'); hold on;
conf = 0.8;
error_ellipse(Sigma1,'mu',mu1, 'conf', conf, 'style', 'k');axis('equal');
error_ellipse(Sigma2,'mu',mu2, 'conf', conf, 'style', 'k');axis('equal');
error_ellipse(Sigma3,'mu',mu3, 'conf', conf, 'style', 'k');axis('equal');
error_ellipse(Sigmak(:,:,1),'mu',mu(:,1), 'conf', conf, 'style', 'r');axis('equal');
error_ellipse(Sigmak(:,:,2),'mu',mu(:,2), 'conf', conf, 'style', 'r');axis('equal');
error_ellipse(Sigmak(:,:,3),'mu',mu(:,3), 'conf', conf, 'style','r');axis('equal');

% EM algorithm - main loop
for i = 1 : NofIter
    
    % E-step
    for k = 1 : K
        Sigmak_inv(:,:,k) = inv(Sigmak(:,:,k));
    end
    
    for k = 1 : K
        for n = 1 : N
            gammakn(k,n)= Pk(k) * ((2 * pi)^(- .5 * l) * det(Sigmak(:,:,k))^( - .5) * exp(-.5 * (X(:,n) - mu(:,k)).' * Sigmak_inv(:,:,k) * (X(:,n) - mu(:,k))));
        end
    end
    gammakn = gammakn ./ repmat(sum(gammakn,1),[K 1]);
    
    % M-step
    Nk = sum(gammakn,2);
    
    for k = 1 : K
        if Nk(k)
            mu(:,k) = 1/Nk(k) * (X * gammakn(k,:).');
            tmp = 1e-5 * eye(l);
            for n = 1 : N
                tmp = tmp + gammakn(k,n) * ((X(:,n) - mu(:,k)) * (X(:,n) - mu(:,k)).');
            end
            Sigmak(:,:,k) = 1/Nk(k) * tmp;
        end
    end
    
    Pk = Nk / N;
    
    for k = 1 : K
        Sigmak_inv(:,:,k) = inv(Sigmak(:,:,k));
    end
    
    % compute the log likelihood
    tmp1 = zeros(N,1);
    for n = 1 : N
        tmp2 = zeros(K,1);
        for k = 1 : K
            tmp2(k) = Pk(k)*(1/(2*pi) * 1/sqrt(det(Sigmak(:,:,k)))*exp(-.5 * (X(:,n) - mu(:,k)).' * Sigmak_inv(:,:,k) * (X(:,n) - mu(:,k))));
        end
        tmp1(n) = log(sum(tmp2));
    end
    lg_lklhd(i) = sum(tmp1);
    
    % plot estimations at iterations 5 and 30
    if i == 5 || i == 50
        figure;
        plot(x1(:,1),x1(:,2),'.k',x2(:,1),x2(:,2),'+k',x3(:,1),x3(:,2),'*k'); hold on;
        error_ellipse(Sigma1,'mu',mu1, 'conf', conf, 'style', 'k');axis('equal');
        error_ellipse(Sigma2,'mu',mu2, 'conf', conf, 'style', 'k');axis('equal');
        error_ellipse(Sigma3,'mu',mu3, 'conf', conf, 'style', 'k');axis('equal');
        error_ellipse(Sigmak(:,:,1),'mu',mu(:,1), 'conf', conf, 'style', 'r');axis('equal');
        error_ellipse(Sigmak(:,:,2),'mu',mu(:,2), 'conf', conf, 'style', 'r');axis('equal');
        error_ellipse(Sigmak(:,:,3),'mu',mu(:,3), 'conf', conf, 'style','r');axis('equal');
    end
end

% plot log likelihood vs iterations
figure; plot(lg_lklhd); ylabel('Likelihood'); xlabel('Iterations');