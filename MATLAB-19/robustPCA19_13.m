%---------------------------------------------------------------------------------------------------------------------------------------------------
% Exercise 19_13
% Robuast PCA.
% Download escalator data from http://perception.i2r.a-star.edu.sg/bk_model/bk_index.htmlhttp://perception.i2r.a-star.edu.sg/bk_model/bk_index.html
%----------------------------------------------------------------------------------------------------------------------------------------------------

load escalator_data 
X = double(X);
nFrames = size(X,2);

% Download approximate proximal gradient method from http://perception.csl.illinois.edu/matrix-rank/sample_code.html
% and add it in the matlab path
addpath('apg')
addpath('apg_partial')
addpath('apg_partial\propack')


maxIter =100;
lambda = 0.01;
%[LL,SS,numIter] = proximal_gradient_rpca(X, lambda,maxIter);

[LL,SS,numIter] = partial_proximal_gradient_rpca(X, lambda,maxIter);


mat  = @(x) reshape( x, m, n );
figure();
colormap( 'Gray' );
k = 1;
for k = 1:nFrames
    
    imagesc([ mat(X(:,k)), mat(LL(:,k)),    mat(SS(:,k)) ] );
    axis off
    axis image
    
    drawnow;
    pause(.05); 
    
end

