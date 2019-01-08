%-----------------------------------------------------------------
%  Exercise 4.25
%  Image Debluring
%  Use MATLAB's Wiener function.
%-----------------------------------------------------------------




clear;
rseed	= 1;
rand('state',rseed);

I =im2double( imread('boat_original.gif') );
[M, N] = size(I);

%add motion using a Point-spread function
PSF = fspecial('motion',20,45);
J = imfilter(I,PSF,'conv', 'circ');


%add noise
noise_mean=0;
noise_var = 0.000001;
J1 = imnoise(J, 'gaussian', noise_mean, noise_var);

%compute root mean square error and psnr
rmse1 = sqrt( sum(sum((I-J1).^2))/M/N );
psnr1 = 20*log10(1/rmse1);


%Wiener filtering (deblurring)
%you can experiemnt with the value of nsr 
nsr = 10^(-2);  %10^-4,   10^-5   10^-6   10^-3   10^-2
K = deconvwnr(J1, PSF, nsr);


%compute root mean square error and psnr
%between original and deblurred pictures
rmse2 = sqrt( sum(sum((I-K).^2))/M/N );
psnr2 = 20*log10(1/rmse2);


%plot results
figure(1);
subplot(2,2,1);
imshow(I);

subplot(2,2,2);
imshow(J1);
title(['PSNR = ', num2str(psnr1,4), 'dB']);

subplot(2,2,3);
imshow(K);
title(['PSNR = ', num2str(psnr2,4), 'dB']);