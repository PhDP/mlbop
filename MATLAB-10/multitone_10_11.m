%--------------------------------------------------------------------------
% Exercise 10.11
% There is a typo in the text concerning this exercise. The parenthesis 
% should read: "( use the “idct.m” Matlab function to compute the inverse 
% DCT transform).
%-------------------------------------------------------------------------


function Matlab10_11

% If you want to use sparselasso.m, SparseLab (https://sparselab.stanford.edu/) need to be downloaded and
% added in matlab path.
% addpath('SparseLab2.1-Core\solvers')

% However, to facilitate the demonstration of this test, a sparse coding algorithm (the OMP that you are asked to
% develope in excercise 10_12) is included at the end of this .m file.


% signal length
l = 2^8;

% number of nonzero frequency components
k = 3;

% number of observations to make
N = 30;


posK=[4,10,30];  

a=[0.3,1,0.75];

n=0:l-1;
theta=zeros(l,1);

% Construct the multitone signal
for i=1:k
   theta=theta+a(i)*cos((pi*(2*posK(i)-1)*n')/(2*l));
end


% question (a)


% The specific signal is sparse in the idct domain (it has a sparse inverse
% dct transform). Here, the original dct.m Matlab file is assumed to be used. 
% Be aware that there are several types of DCT, so in the case that
% you have already include in the matlab path any package/toolbox which have its
% own DCT then the results might not be the expected.
X = idct(theta);

% Instead of using idct, you can work with the DCT matrix as follows
% Phi = dctmtx(l);
% X = Phi'*theta;

figure;subplot(1,2,1);plot(theta);title('Signal');subplot(1,2,2);stem(X);title('DCT domain')


% Question (b)
% Construct the sensing matrix with variance N(0,1/N)
A = randn(N,l)*sqrt(1/N); 
y = A*theta;

% Since it is sparse in the IDCT domain, i.e. A*theta = A*Phi*X = AF*X,
% where X sparse,  AF = A*Phi and Phi is the DCT matrix, Phi = dctmtx(l);.
% Equivalently, using idct (for faster computation than with the DCT matrix), AF is computed as:
AF = idct(A')';


% Use solvelasso.m
% solsA = SolveLasso(AF, y);
% or
% use OMP.m with the sparsity level, i.e. 3, given as input.
solsA = OMP(AF,y,3);
%Take the IDCT (i.e. the DCT) in order to compute the estimated signal. 
theta_hat = dct(solsA);

figure;subplot(1,2,1);plot(theta);title('Original'); subplot(1,2,2);plot(theta_hat);title('Estimated (random sensing matrix)'); 


% Construct the sensing matrix of question (c)
%Question (c)
positions = randperm(l,N);
B = zeros(N,l);
for i=1:N
B(i,positions(i))=1;
end
y = B*theta;


figure; plot(theta); 
hold on
plot(positions,y,'r.')
title('Samples taken')

% Since it is sparse in the IDCT domain, i.e. B*theta = B*Phi*X = BF*X,
% where X sparse,  BF = B*Phi; and Phi is the DCT matrix, Phi = dctmtx(l);.
% Equivalently, using idct (for faster computation than with the DCT matrix), AF is computed as:
BF = idct(B')';

% Use solvelasso.m
% solsB = SolveLasso(BF, y);
% or
% use OMP.m with the sparsity level, i.e. 3, given as input.
solsB = OMP(BF,y,3);

%Take the inverse IDCT (i.e. the DCT) in order to compute the estimated signal. 
theta_hat = dct(solsB);

figure;subplot(1,2,1);plot(theta);title('Original');subplot(1,2,2);plot(theta_hat);title('Estimated (using randomly picked samples)');

end

function theta=OMP(X,y,k)

    residual=y;
    S=zeros(k,1);
    normx = sqrt(sum(X.^2,1))';
    for i=1:1:k,
        proj=X'*residual;
        proj = proj./normx;
        [~,pos]=max(abs(proj));
        pos=pos(1);
        S(i)=pos;
        theta_=pinv(X(:,S(1:i)))*y;

        theta = zeros(size(X,2),1);
        theta(S(1:i)) = theta_;
        residual=y-X*theta;
        
      %  residual=y-X(:,indx(1:i))*theta_; %a faster implementation
    end;
   % theta = zeros(size(X,2),1);
% theta(indx) = theta_;
end