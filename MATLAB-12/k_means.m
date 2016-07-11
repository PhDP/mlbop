function [theta,bel,J]=k_means(X,theta) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%  [theta,bel,J]=k_means(X,theta)
% This function implements the k-means algorithm, which requires
% as input the number of clusters underlying the data set. The algorithm
% starts with an initial estimation of the cluster representatives and
% iteratively tries to move them into regions that are "dense" in data
% vectors, so that a suitable cost function is minimized. This is
% achieved by performing (usually) a few passes on the data set. The
% algorithm terminates when the values of the cluster representatives
% remain unaltered between two successive iterations.
%
% INPUT ARGUMENTS:
%  X:       lxN matrix, each column of which corresponds to
%           an l-dimensional data vector.
%  theta:   a matrix, whose columns contain the l-dimensional (mean)
%           representatives of the clusters.
%
% OUTPUT ARGUMENTS:
%  theta:   a matrix, whose columns contain the final estimations of
%           the representatives of the clusters.
%  bel:     N-dimensional vector, whose i-th element contains the
%           cluster label for the i-th data vector.
%  J:       the value of the cost function (sum of squared Euclidean
%           distances of each data vector from its closest parameter
%           vector) that corresponds to the  estimated clustering.
%
% (c) 2010 S. Theodoridis, A. Pikrakis, K. Koutroumbas, D. Cavouras
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[l,N]=size(X);
[l,m]=size(theta);
e=1;
iter=0;
while(e~=0)
    iter=iter+1;
    theta_old=theta;
    dist_all=[];
    for j=1:m
        dist=sum(((ones(N,1)*theta(:,j)'-X').^2)');
        dist_all=[dist_all; dist];
    end
    [q1,bel]=min(dist_all);
    J=sum(min(dist_all));
    
    for j=1:m
        if(sum(bel==j)~=0)
            theta(:,j)=sum(X'.*((bel==j)'*ones(1,l))) / sum(bel==j);
        end
    end
    e=sum(sum(abs(theta-theta_old)));
end
