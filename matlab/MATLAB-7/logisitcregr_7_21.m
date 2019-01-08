%-----------------------------------------------------------------
%  MATLAB code for Exercise 7.21
%  Logistic Regression
%-----------------------------------------------------------------


clear
format compact
close all

rand('seed',0)

randn('seed',0)
% Definition of mu's and Sigma
% Mean vectors and covariance matrix
m1=[0 2]';  m2=[0 0]'; S1=[4 1.8; 1.8 1]; S2= [4 -1.8; -1.8 1];

% Number of data points
n_points_per_class=1500;

% (i) Data point generation
% Training set
X=[mvnrnd(m1',S1,n_points_per_class); mvnrnd(m2',S2,n_points_per_class)]';
y=[0*ones(1,n_points_per_class) 1*ones(1,n_points_per_class)];
[l,p]=size(X);
%Plot the data set
figure; plot(X(1,y==0),X(2,y==0),'.b',X(1,y==1),X(2,y==1),'.r'); axis equal

% Test set
X_test=[mvnrnd(m1',S1,n_points_per_class); mvnrnd(m2',S2,n_points_per_class)]';
y_test=[0*ones(1,n_points_per_class) 1*ones(1,n_points_per_class)];
[l,p_test]=size(X_test);

%Plot the data set
figure; plot(X_test(1,y_test==0),X_test(2,y_test==0),'.b',X_test(1,y_test==1),X_test(2,y_test==1),'.r'); axis equal


% (ii) Bayes classification of X_test
% Estimation of a priori probabilities
P1=n_points_per_class/p;
P2=P1;
% Estimation of pdf's for each data point
for i=1:p
    p1(i)=(1/(2*pi*sqrt(det(S1))))*exp(-(X_test(:,i)-m1)'*inv(S1)*(X_test(:,i)-m1));
    p2(i)=(1/(2*pi*sqrt(det(S2))))*exp(-(X_test(:,i)-m2)'*inv(S2)*(X_test(:,i)-m2));
end
% Classification of the data points
for i=1:p
    if(P1*p1(i)>P2*p2(i))
        class_test(i)=0;
    else
        class_test(i)=1;
    end
end

% Error probability estimation
Pe=0; %Probability of error
for i=1:p
    if(class_test(i)~=y_test(i))
        Pe=Pe+1;
    end
end
Pe=Pe/p


% (iii) Applying logistic regression, using a gradient descent scheme
X=[X; ones(1,p)];
[l,p]=size(X);

X_test=[X_test; ones(1,p)];
[l,p_test]=size(X_test);

rho=0.001;

theta_ini=[-2 -1 1]';  %rand(l,1)
theta=theta_ini;
e_thres=10e-4;  % threshold for the termination of the algorithm
iter=0;   % Iteration counter
e=1; % The difference between the previous and the current estimate of theta
while e>e_thres
    iter=iter+1; %Increment the iteration counter
    theta_old=theta;  %Store the current theta
    
    s=1./(1+exp(-theta'*X)); %Computation of the vector 's'
    theta=theta-rho*X*(s-y)'; %Updating of the vector 'theta'
    e=sum(abs(theta-theta_old)); %Absolute difference between the current and the previous values of 'theta'
end

% Evaluating the performance of the logistic regression
s_test=1./(1+exp(-theta'*X_test));
Pe_log=0;
for i=1:p_test
    if( ((s_test(i)>0.5)&(y_test(i)~=1)) | ((s_test(i)<0.5)&(y_test(i)~=0)) )
        Pe_log=Pe_log+1;
    end
end
Pe_log=Pe_log/p_test