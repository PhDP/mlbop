%----------------------------------------------------------------------------------------------------------
% Exercise 19.12
% Dictionary learning.
% K-SVD is assumed to be downloaded and properly installed/compiled. The
% version with improvements found here: http://www.cs.technion.ac.il/~ronrubin/Software/ksvdsbox11.zip
% is prefered.
%----------------------------------------------------------------------------------------------------------



clear; clc;



load boats
x_exact = boats;
sigma = 20;
x_noise = x_exact + randn(size(x_exact))*sigma;


blocksize=12; % block size
K=14^2; % number of atoms in the dictionary


X=im2col(x_noise,[blocksize,blocksize],'sliding');

% Training of the  dictionary

params.data = X;
%params.Edata = 0.1; %
params.Tdata = 5;

params.dictsize = K;
params.iternum = 100;
params.memusage = 'high';
params.initdict = Dict_fixed;

[Dict_ksvd,g,err] = ksvd(params,'');

% denoising with ksvd dictionary

X_reconst = zeros(size(X));
subproblems = size(X,2);
 h = waitbar(0,'Please wait...');
for i=1:size(X,2)
    if mod(i,100)==0
    waitbar(i/subproblems,h)
    end

G = omp(Dict_ksvd,X(:,i),[],5);

X_reconst(:,i) = Dict_ksvd*G;
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
recovered_boat_ksvd=yout./Weight; 

figure
subplot(1,3,1)
imshow(boats,[0 255])
x_exact=boats;
subplot(1,3,2)
imshow(x_noise,[0 255])
subplot(1,3,3)
imshow(recovered_boat_ksvd,[0 255])


20*log10(255 * sqrt(numel(x_exact)) / norm(x_exact(:)-recovered_boat_ksvd(:)))

save

