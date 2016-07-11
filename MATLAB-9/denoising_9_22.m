%--------------------------------------------------------------------
% Exercise 9.22
% SparseLab is used here (downloaded from https://sparselab.stanford.edu/)
% It can be replaced with any other sparse construction tool/algorithm.
%------------------------------------------------------------------------

function Matlab9_22



addpath(genpath('C:\research\pocs\SparseLab2.1-Core'));


load boats
x_exact = boats;
sigma = 20;
x_noise = x_exact + randn(size(x_exact))*sigma;


blocksize=12; % block size
K=14^2; % number of atoms in the dictionary


X=im2col(x_noise,[blocksize,blocksize],'sliding');

%Question (a)
% use fixed dictionary
n = 0:blocksize-1;
DCT=zeros(blocksize,sqrt(K));
for k=0:1:sqrt(K)-1,
   V=cos(n'*k*pi/sqrt(K));
   
   DCT(:,k+1)=V/norm(V);
end;

% Question (b)
DCT=kron(DCT,DCT);
Dict_fixed = DCT;
for k=1:size(Dict_fixed,2)
   Dict_fixed(:,k) = Dict_fixed(:,k)/norm(Dict_fixed(:,k)); 
end


% denoising with DCT

X_reconst = zeros(size(X));
subproblems = size(X,2);
 h = waitbar(0,'Please wait...');
for i=1:size(X,2)
    if mod(i,100)==0
    waitbar(i/subproblems,h)
    end
    
    G = SolveLasso(Dict_fixed, X(:,i),size(Dict_fixed,2),'lars',20+1);
% It will be faster using OMP (check
% http://www.cs.technion.ac.il/~ronrubin/Software/ompbox10.zip for a fast
% implementation)
%G = omp(Dict_fixed,X(:,i),[],5);
% or use the OMP you will write in Matlab excersize 10.12
% G = OMP(Dict_fixed,X(:,i),5);

X_reconst(:,i) = Dict_fixed*G;
end
close(h)


% ========================================================
%this is contained in Chapter_15_CoreInpaining1.m of Elad's book
%"Sparse and Redundant Representations: From Theory to Applications in Signal and Image Processing" , Springer, 2010
% ========================================================
N=size(x_exact,1); 
n=blocksize; %sqrt(size(D,1)); 
yout=zeros(N,N); 
Weight=zeros(N,N); 
i=1; j=1;
for k=1:1:(N-n+1)^2,
    patch=reshape(X_reconst(:,k),[n,n]); 
    yout(i:i+n-1,j:j+n-1)=yout(i:i+n-1,j:j+n-1)+patch; 
    Weight(i:i+n-1,j:j+n-1)=Weight(i:i+n-1,j:j+n-1)+1; 
    if i<N-n+1 
        i=i+1; 
    else
        i=1; j=j+1; 
    end;
end;
recovered_boat_dct=yout./Weight; 


20*log10(255 * sqrt(numel(x_exact)) / norm(x_exact(:)-x_noise(:)))
20*log10(255 * sqrt(numel(x_exact)) / norm(x_exact(:)-recovered_boat_dct(:)))

figure
subplot(1,3,1)
imshow(boats,[0 255])
x_exact=boats;
subplot(1,3,2)
imshow(x_noise,[0 255])
subplot(1,3,3)
imshow(recovered_boat_dct,[0 255])

end

function theta=OMP(X,y,S)

    residual=y;
    indx=zeros(S,1);
    normx = sqrt(sum(X.^2,1))';
    for i=1:1:S,
        proj=X'*residual;
        proj = proj./normx;
        [~,pos]=max(abs(proj));
        pos=pos(1);
        indx(i)=pos;
        theta_=pinv(X(:,indx(1:i)))*y;
        theta = zeros(size(X,2),1);
        theta(indx(1:i)) = theta_;

        residual=y-X*theta;
      %  residual=y-X(:,indx(1:i))*theta_;
    end;
   % theta = zeros(size(X,2),1);
% theta(indx) = theta_;
end

