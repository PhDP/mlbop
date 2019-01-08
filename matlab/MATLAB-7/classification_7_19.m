%-----------------------------------------------------------------
%  MATLAB code for Exercise 7.19
%  Bayesian Classification
%-----------------------------------------------------------------


clear
format compact
close all

randn('seed',0)
% Definition of mu's and Sigma
% Mean vectors and covariance matrix
m1=[0 0]';  m2=[2 2]'; S=[1 .25; .25 1];
% Number of data points
n_points_per_class=500;

% (i) Data point generation
X=[mvnrnd(m1',S,n_points_per_class); mvnrnd(m2',S,n_points_per_class)]';
label=[ones(1,n_points_per_class) 2*ones(1,n_points_per_class)];
[l,p]=size(X);
%Plot the data set
figure; plot(X(1,label==1),X(2,label==1),'.b',X(1,label==2),X(2,label==2),'.r')


% (ii) Bayes classification of X
% Estimation of a priori probabilities
P1=n_points_per_class/p;
P2=P1;
% Estimation of pdf's for each data point
for i=1:p
    p1(i)=(1/(2*pi*sqrt(det(S))))*exp(-(X(:,i)-m1)'*inv(S)*(X(:,i)-m1));
    p2(i)=(1/(2*pi*sqrt(det(S))))*exp(-(X(:,i)-m2)'*inv(S)*(X(:,i)-m2));
end
% Classification of the data points
for i=1:p
    if(P1*p1(i)>P2*p2(i))
        class(i)=1;
    else
        class(i)=2;
    end
end

% (iii) Error probability estimation
Pe=0; %Probability of error
for i=1:p
    if(class(i)~=label(i))
        Pe=Pe+1;
    end
end
Pe=Pe/p

%Plot the data set
figure; plot(X(1,class==1),X(2,class==1),'.b',X(1,class==2),X(2,class==2),'.r')


% Definition of the loss matrix
L=[0 1; .005 0];

% (iv) Classification of data points according to the average risk
% minimization rule
% Estimation of pdf's for each data point
for i=1:p
    p1(i)=(1/(2*pi*sqrt(det(S))))*exp(-(X(:,i)-m1)'*inv(S)*(X(:,i)-m1));
    p2(i)=(1/(2*pi*sqrt(det(S))))*exp(-(X(:,i)-m2)'*inv(S)*(X(:,i)-m2));
end
% Classification of the data points
for i=1:p
    if(L(1,2)*P1*p1(i)>L(2,1)*P2*p2(i))
        class_loss(i)=1;
    else
        class_loss(i)=2;
    end
end

% (v) Error probability estimation
Ar=0;  %Average risk
for i=1:p
    if(class_loss(i)~=label(i))
        if(label(i)==1)
            Ar=Ar+L(1,2);
        else
            Ar=Ar+L(2,1);
        end
    end
end
Ar=Ar/p

%Plot the data set
figure; plot(X(1,class_loss==1),X(2,class_loss==1),'.b',X(1,class_loss==2),X(2,class_loss==2),'.r')


% The average risk minimization rule leads to smaller value of the average risk compared to the 
% probability of error obtained by the classic Bayes classification rule.
% In the former case the classification rule decides for almost all the
% overapping region between the two classes, in favor of \omega_1. This is
% due to the fact that a classification error on data steming from \omega_2 is ``cheap'',
% compared to an opposite classification error.