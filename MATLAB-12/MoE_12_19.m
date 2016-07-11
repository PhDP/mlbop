%-----------------------------------------------------------------
%  Exercise 12.19
%  Mixture of experts
%-----------------------------------------------------------------


clear; close all; clc; format long eng; format compact;

rng('default');

% max iterations
Imax = 100; 

% # of models
K = 2; 

% # of input points
N = 50; 

% Generate data
x = linspace(-1,1,N);
idx(x < -0.5) = 1;
idx( (x > -0.5) & (x < 0.5) ) = 2;
idx(x > 0.5) = 1;
idx1 = (idx==1);
idx2 = (idx==2);
w = 0.01*randn(1,K);
b = [-1 1];
t = x.*w(idx) + b(idx) + 0.1*randn(1,N);

x = [x.', ones(N,1)];

% EM algorithm
% initialization
w_em = randn(2,K) - rand(2,K);
gnk = 1/K * ones(N,K) + .02 * randn(N,K); 
pk = 1/K * ones(K,1) ;  
vita = 10;
i = 0;
for i = 1 : Imax
        
    % plot the results at specified iterations
    if i == 1 || i == 4 || i == 100
        figure; hold on; plot(x(idx1,1),t(idx1),'ro'); plot(x(idx2,1),t(idx2),'ko'); plot(x(:,1),x(:,1).*repmat(w_em(1,1),N,1) + repmat(w_em(2,1),N,1),'r-');plot(x(:,1),x(:,1).*repmat(w_em(1,2),N,1) + repmat(w_em(2,2),N,1),'k-'); box on; hold off; 
        figure; line([x(:,1).'; x(:,1).'], [zeros(1,N); ones(1,N)],'Color','k'); hold on; stem(x(:,1),gnk(:,1),'rx'); box on; hold off;
    end
    
    % E-step
    for n = 1 : N
        tmp = 0;
        for j = 1 : K
            tmp = tmp + pk(j) * sqrt(.5 / pi) * sqrt(vita) * exp(-.5 * vita * (t(n) - w_em(:,j).' * x(n,:).')^2);
        end
        
        for k = 1 : K
            gnk(n,k) = pk(k) * sqrt(.5 / pi) * sqrt(vita) * exp(-.5 * vita * (t(n) - w_em(:,k).' * x(n,:).')^2) / tmp;
        end
    end
    
    pk = 1/N * sum(gnk,1);
    
    % M-step
    w_em_old = w_em;
    for k = 1 : K
        w_em(:,k) = (x.' * diag(gnk(:,k)) * x) \ (x.' * diag(gnk(:,k)) * t.');
    end
    
    tmp = 0;
    for n = 1 : N
        for k = 1 : K
            tmp = tmp + gnk(n,k) * (t(n) - w_em(:,k).' * x(n,:).')^2;
        end
    end
    vita = N / tmp;
    
end

