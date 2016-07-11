%-----------------------------------------------------------------
%  Exercise 12.18
%  Gaussian mixture models: k-means and EM comparison
%  Use kmeans function
%-----------------------------------------------------------------


clc; clear; close all; format compact; format long eng;

% data points
N1 = 100; 
N2 = 20; % or (100)
N = N1 + N2;

% #classes (try different numbers)
K = 2; 

% data dimension
l = 2; 

rng('default'); 

% generate samples out of two Gaussians
mu = [ .9 -1.2;1.02 -1.3];
Sigma1 = [.5 .081;.081 .7]; 
Sigma2 = [.4 .02;.02 .3]; 
x1 = mvnrnd(mu(:,1), Sigma1, N1);
x2 = mvnrnd(mu(:,2), Sigma2, N2);

% set marker color
marker_color = [1 0 0; 0 0 0; .5 .5 .5];

% plot generated data
figure; axis equal; box on; 
hold on;
plot(x1(:,1),x1(:,2),'.','MarkerSize',10','Color',marker_color(2,:));
plot(x2(:,1),x2(:,2),'.','MarkerSize',10','Color',marker_color(1,:));
hold off;

% data matrix
X = [x1.' x2.'];

% kmeans algorithm
theta = zeros(l,K);
[theta, bel, L] = k_means(X,theta);

% plot kmeans clustering
figure; axis equal; box on; 
hold on;
for k = 1 : K
    plot(X(1,bel==k),X(2,bel==k),'.','MarkerSize',10','Color',marker_color(k,:))
end
hold off;

% EM algorithm
% initialization
NofIter = 300;
conf = .8;
Pk = 1/K* ones(K,1); 
mu = randn(l,K);
Sigmak = repmat( rand(2) .* eye(l,l), [1 1 K]);
Sigmak_inv = Sigmak;
gammakn = zeros(K,N);

% plot initial estimates
figure; axis equal; box on; 
hold on;
plot(x1(:,1),x1(:,2),'.','MarkerSize',10','Color',marker_color(2,:));
plot(x2(:,1),x2(:,2),'.','MarkerSize',10','Color',marker_color(1,:));
for k = 1 : K,  if Pk(k)~=0, error_ellipse((Sigmak(:,:,k)),'mu',mu(:,k), 'conf', conf,'style','r');  end, end
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
    
    for k = 1 : K
        Sigmak_inv(:,:,k) = inv(Sigmak(:,:,k));
    end
end

% plot final estimates
figure; axis equal; box on; 
hold on;
plot(x1(:,1),x1(:,2),'.','MarkerSize',10','Color',marker_color(2,:));
plot(x2(:,1),x2(:,2),'.','MarkerSize',10','Color',marker_color(1,:));
for k = 1 : K,  if Pk(k)~=0, error_ellipse(Sigmak(:,:,k),'mu',mu(:,k), 'conf', conf,'style','r'); axis('equal'); end, end
hold off;
