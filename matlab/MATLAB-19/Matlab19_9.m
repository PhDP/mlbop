%-------------------------------------------------------------------------------
% Exercise 19.9
% You must download the images from
% http://cgi.di.uoa.gr/~stheodor/faces.rar
%--------------------------------------------------------------------------------

function Matlab19_9

sizeimages = 168;

dirname = [pwd,filesep,'faces'];


img_count = 0;

% Read the face images
files_struct = dir(dirname);
X = zeros(sizeimages^2,length(files_struct)-2);
for i = 3:length(files_struct),
    filename = files_struct(i).name;
    full_pathname = [dirname, '/', filename];
    img = imread(full_pathname);
    img_count = img_count + 1;
    X(:,img_count) = uint8(img(:));
end

N = size(X,2);

%  Subtract out the means from each parameter
mX = mean(X,2);
X = X - repmat(mX,1,N);


% Effectively compute the eigenvectors based on X'X.
[S,eigenvalues,~] = svd(X'*X,0);
eigenvectors = X*S;
for i = 1:size(X,2)
    eigenvectors(:,i) = 1/sqrt(eigenvalues(i,i))*eigenvectors(:,i);
end

% Plot some eigenfaces. Change some signs for fun.
sp = [-1 -2 -6 -7 8 -10 -11 17];
pp =0;
figure
for ii=1:8
    img = reshape(sign(sp(ii))*eigenvectors(:,abs(sp(ii))), sizeimages, sizeimages);
    pp=pp+1;
    subplot(2,4,pp)
    imagesc(img)
    colormap gray
    axis square
    axis off;
end

% choose a face randomly
pickrandom = randperm(N,1);
testimg = X(:,pickrandom);

subim = 0;
subim=subim+1;
figure;
subplot(1,5,subim)
img = reshape(testimg+mX, 168, 168); %The mean is added back again in order to display it correctly
imagesc(img); colormap gray; axis square; axis off
title('original')

% Reconstruct using different number of eigenvectors
PP = eigenvectors'*(double(testimg));

[a,b] = sort(abs(PP),'descend');

subim=subim+1;
subplot(1,5,subim)
howmany = 5;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square; axis off;
title([num2str(howmany), ' eigen'])

subim=subim+1;
subplot(1,5,subim)
howmany = 30;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square; axis off;
title([num2str(howmany), ' eigen'])

subim=subim+1;
subplot(1,5,subim)
howmany = 100;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square;  axis off;
title([num2str(howmany), ' eigen'])

subim=subim+1;
subplot(1,5,subim)
howmany = 600;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square ;  axis off;
title([num2str(howmany), ' eigen'])

end


