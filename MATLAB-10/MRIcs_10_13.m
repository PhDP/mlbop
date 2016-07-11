%-------------------------------------------------------------
% Exercise 10.13
% MRI reconstruction.
% The NESTA toolbox need to be downloaded and added to the matlab path 
% (http://statweb.stanford.edu/~candes/nesta/nesta.html).
% Hadamard function (distributed with NESTA) might also need compiling.
%---------------------------------------------------------------

function Matlab10_13

        
        
        addpath(genpath('NESTA_v1.1'))
        


% Load The phantom
load phant
X = ph;

x = X(:);
n = 256;
% number of radial lines in the Fourier domain
Ls = 17;

% Take radial lines 
[~,~,~,OMEGA] = LineMask(Ls,n);


% Perform NESTA (this part resemples the DemoTV.m demonstration file in NESTA 
A = @(z) A_fhp(z, OMEGA);
At = @(z) At_fhp(z, OMEGA, n);

b = A(x);

Psi = @(x) x; PsiT = Psi;

        
U = PsiT;
Ut = Psi;

mu = 0.0005; %--- can be chosen to be small
opts = [];
opts.maxintiter = 5;
opts.TOlVar = 1e-5;
opts.verbose = 1;
opts.maxiter = 5000;
opts.U = U;
opts.Ut = Ut;
opts.stoptest = 1;  
opts.typemin = 'tv';
delta = 0.001; 

[x_nesta,niter,resid,err] = NESTA(A,At,b,mu,delta,opts);
Xnesta = reshape(x_nesta,n,n);
figure; imshow(Xnesta)
 
backproj = At(b);
backproj = reshape(backproj,n,n);
figure;imshow(backproj)
   


end




