%-----------------------------------------------------------------
%  Exercise 13.27
%  Hyperspectral unmixing
%  Notice that the unmixing process, using the full image,
%  requires a few hours. To shorten the time, part of the image 
%  can be used.
%-----------------------------------------------------------------


clc ; clear; close all; format compact; format long eng;

% load the hyperspectral image
load cuprite_ref.mat

% number of endmembers
Nen = 14; 

%  number of pixels - to shorten computation time, select 100 pixels only
Np = size(x,2); 

% extract endmembers using VCA algorithm
Phi = VCA(x,'Endmembers',Nen);
Phi_gram = Phi.' * Phi;

[M,N] = size(Phi);

const = 1e-10;
MaxIter = 10000;

w_vb = zeros(N,Np);

for n = 1 : Np
    % variational Bayes abundance estimation 
    w_vb(:,n) = mpLaplace_truncatedGaussian_vB(Phi, x(:,n), Phi_gram, MaxIter, const);
    
end

% reshape variables to display abundance maps
im = zeros(Lines, Columns, L); 
w_vb_im = zeros(Lines, Columns, Nen); 

for i = 1 : Lines, 
    for j = 1 : Columns, 
        im(i,j,:) = x(:,(j-1) * Lines + i ); 
        w_vb_im(i,j,:) = w_vb(:,(j-1) * Lines + i ); 
    end,
end

% display abundance maps
for i = 1 : Nen, figure; imagesc(w_vb_im(:,:,i),[0 1]); axis('image'); axis('off'); end

