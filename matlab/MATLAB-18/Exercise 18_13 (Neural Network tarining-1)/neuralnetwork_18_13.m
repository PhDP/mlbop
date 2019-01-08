%-----------------------------------------------------------------
%  Exercise 18.13
%  Feedforward Neural Networks
%-----------------------------------------------------------------

clear
format compact
close all

randn('seed',0);
% 1. Generate data set X1 
l=2; %Dimensionality
m1=[-5 5; 5 -5]'; % centroids
m2=[-5 -5; 0 0; 5 5]';
[l,c1]=size(m1); %no of gaussians per class
[l,c2]=size(m2);

P1=ones(1,c1)/c1; % weights of the mixture model
P2=ones(1,c2)/c2;
s=1; % variance

% Generate the training data from the first class
N1=100; %60; %Number of first class data points
for i=1:c1
    S1(:,:,i)=s*eye(l);
end
sed=0; %Random generator seed
[class1_X,class1_y]=mixt_model(m1,S1,P1,N1,sed);

% Generate the training data from the second class
N2=150; %80; %Number of second class data points
for i=1:c2
    S2(:,:,i)=s*eye(l);
end
sed=0; 
[class2_X,class2_y]=mixt_model(m2,S2,P2,N2,sed);

% Form X1
X1=[class1_X  class2_X]; %Data vectors
y1=[ones(1,N1) -ones(1,N2)]; %Class labels
figure(1), hold on
figure(1), plot(X1(1,y1==1),X1(2,y1==1),'r.',X1(1,y1==-1),X1(2,y1==-1),'bx')

% Generate test set X2

% Data of the first class
sed=100; %Random generator seed. This time we set this value to 100
[class1_X,class1_y]=mixt_model(m1,S1,P1,N1,sed);

% Data of the second class
sed=100; %Random generator seed
[class2_X,class2_y]=mixt_model(m2,S2,P2,N2,sed);

%Production of the unified data set
X2=[class1_X class2_X]; %Data vectors
y2=[ones(1,N1) -ones(1,N2)]; %Class labels


%Definition and training of the network
rand('seed',100) 
randn('seed',100)

iter=9000; %Number of iterations
code=1; %Code for the chosen training algorithm
k=2; %number of hidden layer nodes
lr=.01; %learning rate
par_vec=[lr 0 0 0 0];
[net,tr]=NN_training(X1,y1,k,code,iter,par_vec);

% Compute the training and the test errors
pe_train=NN_evaluation(net,X1,y1)
pe_test=NN_evaluation(net,X2,y2)

% Plot the data points as well as the decision regions of the FNN
maxi=max(max([X1'; X2']));
mini=min(min([X1'; X2']));
bou=[mini maxi];
fig_num=2; % figure handle
resolu=(bou(2)-bou(1))/100; % figure resolution
plot_NN_reg(net,bou,resolu,fig_num); % plot decision region
figure(fig_num), hold on % plot training set
%figure(fig_num), plot(X1(1,y1==1),X1(2,y1==1),'r.', X1(1,y1==-1),X1(2,y1==-1),'bx')
figure(fig_num), plot(X2(1,y2==1),X2(2,y2==1),'r.', X2(1,y2==-1),X2(2,y2==-1),'bx')

% Plot the training error versus the number of iterations
figure(3), plot(tr.perf)
