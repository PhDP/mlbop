%-------------------------------------------------------------------------------
% Exercise 19.10
% You must download the images from
% http://cgi.di.uoa.gr/~stheodor/faces.rar
%--------------------------------------------------------------------------------

function Matlab19_10

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

% choose a face randomly
N = size(X,2);
pickrandom = randperm(N,1);
testimg = X(:,pickrandom);

subim = 0;
subim=subim+1;
figure;
subplot(1,3,subim)
img = reshape(testimg, 168, 168); 
imagesc(img); colormap gray; axis square; axis off
title('original')

% exclude the selected face from the face-set
X(:,pickrandom) = [];

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



% Reconstruct using different number of eigenvectors
PP = eigenvectors'*(double(testimg));

[a,b] = sort(abs(PP),'descend');

subim=subim+1;
subplot(1,3,subim)
howmany = 300;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square;  axis off;
title([num2str(howmany), ' eigen'])

subim=subim+1;
subplot(1,3,subim)
howmany = 1000;
recimg = eigenvectors(:,b(1:howmany))*PP(b(1:howmany))+mX;
img = reshape(recimg, 168, 168); imagesc(img); colormap gray; axis square ;  axis off;
title([num2str(howmany), ' eigen'])

end


