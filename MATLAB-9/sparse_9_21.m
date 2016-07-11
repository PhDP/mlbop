%--------------------------------------------------------------------------
% Exercise 9.21
% SparseLab is used here (downloaded from https://sparselab.stanford.edu/)
% It can be replaced with any other sparse construction tool/algorithm.
%-------------------------------------------------------------------------

% addpath(genpath('C:\research\pocs\SparseLab2.1-Core'));

N=30;
K=5;
l=100;
rep = 500;

theta = zeros(l,1);
theta(1:K) = randn(K,1);
X = randn(N,l)*(1/sqrt(N));
y = X*theta;

% Question (a)
sols = SolveLasso(X, y);
error =norm(sols-theta);
disp(error)


% Question (b)
Err = zeros(rep,1);
for epan=1:rep

X = randn(N,l)*(1/sqrt(N)); %X=X*diag(1./sqrt(sum(X.^2)));
y = X*theta;

sols = SolveLasso(X, y,length(theta),'lasso');
errorX =norm(sols-theta);
% disp(errorA)
    
Err(epan) = errorX<10^(-8) + 0;
end
probrandn = sum(Err)/rep; 
fprintf('Random Sensing Mtx: %2.2f \n',probrandn)

% Question (c)
Err = zeros(rep,1);
for epan=1:rep

% Construct DCT based sensing matrix
X = dctmtx(l); X=X(randperm(l,N),:);

y = X*theta;

sols = SolveLasso(X, y,length(theta),'lasso');
errorX =norm(sols-theta);
% disp(errorA)
    
Err(epan) = errorX<10^(-8) + 0;
end
probDCT = sum(Err)/rep; 
fprintf('DCT Sensing Mtx: %2.2f \n',probDCT)


% Question (d)
p = [1,9,25,36,64];
Err = zeros(rep,1);
for pval = p
for epan=1:rep

% Construct sparse sensing matrix
OK = false;
while ~OK
    kk = zeros(N*l,1);
    P = randperm(N*l);
    numofzeros = round(N*l*(1-1/sqrt(pval)));
    P(N*l-numofzeros+1:end) = [];
    kk(P) = sqrt( (sqrt(pval)/N) )*ones(length(P),1).*sign(randn(length(P),1));
    X = zeros(N,l);
    X(1:l*N) = kk;
    % Check if it is full rank
    OK = rank(X)==N; 
end

y = X*theta;

sols = SolveLasso(X, y,length(theta),'lasso');
errorX =norm(sols-theta);
% disp(errorA)
    
Err(epan) = errorX<10^(-8) + 0;
end
probsparse = sum(Err)/rep; 
fprintf('Sparse Sensing Mtx, p=%i: %2.2f \n',pval,probsparse)


end

disp('')
