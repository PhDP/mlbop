%-----------------------------------------------------------------
%  MATLAB code for Exercise 7.20
% Naive Bayes Classifier
%-----------------------------------------------------------------


clear
format compact
close all

randn('seed',0)
% Definition of mu's and Sigma
% Mean vectors and covariance matrix
m1=[0 2]';  m2=[0 0]'; S1=[4 0; 0 1]; S2= [4 0; 0 1];  %S1=[4 1.8; 1.8 1]; S2= [4 1.2; 1.2 1];

% Number of data points
n_points_per_class=5000;

% (i) Data point generation
X=[mvnrnd(m1',S1,n_points_per_class); mvnrnd(m2',S2,n_points_per_class)]';
label=[ones(1,n_points_per_class) 2*ones(1,n_points_per_class)];
[l,p]=size(X);
%Plot the data set
figure; plot(X(1,label==1),X(2,label==1),'.b',X(1,label==2),X(2,label==2),'.r'); axis equal


% (ii) Bayes classification of X
% Estimation of a priori probabilities
P1=n_points_per_class/p;
P2=P1;
% Estimation of pdf's for each data point
for i=1:p
    p1(i)=(1/(2*pi*sqrt(det(S1))))*exp(-(X(:,i)-m1)'*inv(S1)*(X(:,i)-m1));
    p2(i)=(1/(2*pi*sqrt(det(S2))))*exp(-(X(:,i)-m2)'*inv(S2)*(X(:,i)-m2));
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
figure; plot(X(1,class==1),X(2,class==1),'.b',X(1,class==2),X(2,class==2),'.r'); axis equal

%(iv) 
% Estimation of pdf's for each data point
for i=1:p
    p1_naive(i)=( (1/(sqrt(2*pi*S1(1,1))))*exp(-(X(1,i)-m1(1))'*(X(1,i)-m1(1))/(2*S1(1,1))) ) * ( (1/(sqrt(2*pi*S1(2,2))))*exp(-(X(2,i)-m1(2))'*(X(2,i)-m1(2))/(2*S1(2,2))) );
    p2_naive(i)=( (1/(sqrt(2*pi*S2(1,1))))*exp(-(X(1,i)-m2(1))'*(X(1,i)-m2(1))/(2*S2(1,1))) ) * ( (1/(sqrt(2*pi*S2(2,2))))*exp(-(X(2,i)-m2(2))'*(X(2,i)-m2(2))/(2*S2(2,2))) );
end
% Classification of the data points
for i=1:p
    if(P1*p1_naive(i)>P2*p2_naive(i))
        class_naive(i)=1;
    else
        class_naive(i)=2;
    end
end

% (v) Error probability estimation
Pe_naive=0; %Probability of error
for i=1:p
    if(class_naive(i)~=label(i))
        Pe_naive=Pe_naive+1;
    end
end
Pe_naive=Pe_naive/p

%Plot the data set
figure; plot(X(1,class_naive==1),X(2,class_naive==1),'.b',X(1,class_naive==2),X(2,class_naive==2),'.r'); axis equal
