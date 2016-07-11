%-----------------------------------------------------------------
%  Exercise 13.22
%  Gaussian mixture modeling with variational EM algorithm
%  Use error_ellipse function.
%-----------------------------------------------------------------


clc; clear; close all; format long eng; format compact;

rng('default');

% #samples
N = 300; 

% dimension
l = 2; 

% #Gaussians
K = 25; 

% true parameters
mu = [-2.5 -4 2 0.1 3;2.5 -2 -1 .2 3];
Sigma1 = [.5 .081;.081 .7]; Sigma2 = [.4 .02;.02 .3]; Sigma3 = [.6 .531;.531 .9]; Sigma4 = [.5 .22;.22 .8]; Sigma5 = [.88 .2;.2 .22];

% generate data
x1 = mvnrnd(mu(:,1), Sigma1, N/5);
x2 = mvnrnd(mu(:,2), Sigma2, N/5);
x3 = mvnrnd(mu(:,3), Sigma3, N/5);
x4 = mvnrnd(mu(:,4), Sigma4, N/5);
x5 = mvnrnd(mu(:,5), Sigma5, N/5);

% plot generated data
% plot(x1(:,1),x1(:,2),'.b',x2(:,1),x2(:,2),'.g',x3(:,1),x3(:,2),'.r',x4(:,1),x4(:,2),'.m',x5(:,1),x5(:,2),'.k'); %plot data
 
% data matrix
X = [x1.' x2.' x3.' x4.' x5.' ];

% EM algorithm

% initialization
NofIter = 300;
conf = .8;
Pk = 1/K* ones(K,1);
mu = 2 * randn(l,K);
Sigmak = rand(l,l,K) + 5e1 *  repmat(eye(l),[1 1 K]);
Sigmak_inv = Sigmak;
gammakn = zeros(K,N);

% plot initial Gaussians
figure; axis('equal');
hold on;
plot(x1(:,1),x1(:,2),'k.',x2(:,1),x2(:,2),'k.',x3(:,1),x3(:,2),'k.',x4(:,1),x4(:,2),'*k',x5(:,1),x5(:,2),'.k'); axis('equal');
for k = 1 : K,  if Pk(k)~=0, error_ellipse((Sigmak(:,:,k)),'mu',mu(:,k), 'conf', conf, 'style', 'r'); axis('equal'); end, end
hold off;

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
end

% plot EM estimates
figure; axis('equal');
hold on;
plot(x1(:,1),x1(:,2),'k.',x2(:,1),x2(:,2),'k.',x3(:,1),x3(:,2),'k.',x4(:,1),x4(:,2),'*k',x5(:,1),x5(:,2),'.k'); axis('equal');
for k = 1 : K,  if Pk(k)~=0, error_ellipse((Sigmak(:,:,k)),'mu',mu(:,k), 'conf', conf, 'style','r'); axis('equal'); end, end 
axis('equal');
hold off;


% Variational Bayes algorithm

% initialization
NofIter = 300;
conf = .8;
beta = 1;
nu = ones(K,1);
W0_inv = eye(l);
EQk = repmat(W0_inv,[1 1 K]);
Pk = 1/K * ones(K,1);
ElnQk = ones(K,1);
mu_tilde =  randn(l,K);
Q_tilde_inv = rand(l,l,K) + 5e1 *  repmat(eye(l),[1 1 K]);
Q_tilde = 1./(Q_tilde_inv);
Emukmuk = Q_tilde;
Wk_inv = zeros(l,l,K);
pik = zeros(K,N);

% plot initial Gaussians
figure; axis('equal');
hold on;
plot(x1(:,1),x1(:,2),'k.',x2(:,1),x2(:,2),'k.',x3(:,1),x3(:,2),'k.',x4(:,1),x4(:,2),'*k',x5(:,1),x5(:,2),'.k'); axis('equal');
for k = 1 : K,  if Pk(k)~=0, error_ellipse((Q_tilde_inv(:,:,k)),'mu',mu_tilde(:,k), 'conf', conf,'style','r'); axis('equal'); end, end
hold off;

for i = 1 : NofIter
    % E-step 1a
    for k = 1 : K
        for n = 1 : N
            pik(k,n) = Pk(k) * exp(.5 * ElnQk(k) - .5 * trace(EQk(:,:,k) * (X(:,n) * X(:,n).' - X(:,n) * mu_tilde(:,k).' - mu_tilde(:,k) * X(:,n).' + Emukmuk(:,:,k)) ) );
        end
    end
    pnk = pik ./ repmat(sum(pik,1), [size(pik,1) 1]);
    
    % E-step 1b
    pnk_sumk = sum(pnk,2);
    for k = 1 : K
        Q_tilde(:,:,k) = beta * eye(l) + EQk(:,:,k) * pnk_sumk(k);
        Q_tilde_inv(:,:,k) = inv(Q_tilde(:,:,k));
        mu_tilde(:,k) =  Q_tilde_inv(:,:,k) * EQk(:,:,k) * (X * pnk(k,:).');
        Emukmuk(:,:,k) = Q_tilde_inv(:,:,k) + mu_tilde(:,k) * mu_tilde(:,k).';
    end
    
    % E-step 1c
    nuk =  nu + pnk_sumk;
    for k = 1 : K
        Wk_inv(:,:,k) = W0_inv;
        for n = 1 : N
            Wk_inv(:,:,k) = Wk_inv(:,:,k) + pnk(k,n) * (X(:,n) * X(:,n).' - mu_tilde(:,k) * X(:,n).' - X(:,n) * mu_tilde(:,k).' + Emukmuk(:,:,k));
        end
        Wk = inv(Wk_inv(:,:,k));
        EQk(:,:,k) = nuk(k) * Wk; 
        ElnQk(k) = psi(.5 * nuk(k)) + psi(.5 * (nuk(k) - 1) ) + l * log(2) + log(det(Wk));
    end
    
    % M-step
    Pk = 1/N * pnk_sumk;
end

% plot the results
figure;
hold on;
plot(x1(:,1),x1(:,2),'k.',x2(:,1),x2(:,2),'k.',x3(:,1),x3(:,2),'k.',x4(:,1),x4(:,2),'*k',x5(:,1),x5(:,2),'.k'); axis('equal');
for k = 1 : K,  if Pk(k)~=0, error_ellipse((inv(EQk(:,:,k))),'mu',mu_tilde(:,k), 'conf', conf,'style','r'); axis('equal'); end, end
axis('equal');
hold off;