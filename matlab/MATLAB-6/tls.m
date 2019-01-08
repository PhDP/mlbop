%=====================================================================
% Exercise 6.22
% This script implements the Total Least Squares algoritm
%=====================================================================


function x = tls(A,b,thresh)

% Solves the linear equation Ax=b using
% truncated total least squares.

[m n] = size(A);
if sum(size(b)-[m 1])
error('A, b size mis-match')
end
if nargin==2
thresh = eps;
end

% augmented matrix
Z = [full(A');b'];
[V W] = svd(Z);% Compute the SVD of the augmented matrix

% find sing val above
d = diag(W);
k = sum(d<thresh);
q = n - k + 1;

 V12 = V(1:n,q:end);
V22 = V(n+1,q:end);
x = -V12 * V22' ./ norm(V22).^2;