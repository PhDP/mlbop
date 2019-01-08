% This is for theodoridis book. Boosting example FOR THE TEXT

clear
format compact
close all

% Generating the training date set
randn('seed',0)
l=20;
m11=[0*ones(1,l/2) 0*ones(1,l/2)]'; m12=[1*ones(1,l/2) 1*ones(1,l/2)]'; m21=[0*ones(1,l/2) 1*ones(1,l/2)]'; %m22=[3*ones(1,l/2) 0*ones(1,l/2)]';
S11=1*eye(l); S12=1*eye(l); S21=1*eye(l); %S22=1*eye(l);
%S11=[.2 0; 0 2]; S12=[3 0; 0 0.5]; S21=[5 0; 0 0.5]; S22=[7 0; 0 0.5]; %S3=[8 0; 0 0.5];
n_of_points_per_group=100;

X=[mvnrnd(m11',S11,n_of_points_per_group); mvnrnd(m12',S12,n_of_points_per_group); mvnrnd(m21',S21,n_of_points_per_group)]'; %mvnrnd(m22',S22,n_of_points_per_group)]'; %mvnrnd(m3',S3,n_of_points_per_group)]';
label=[ones(1,n_of_points_per_group) ones(1,n_of_points_per_group) 2*ones(1,n_of_points_per_group)]; % 2*ones(1,n_of_points_per_group)];  % 3*ones(1,n_of_points_per_group)];
[l,p]=size(X);
%Plot the training data set
%figure; plot(X(1,label==1),X(2,label==1),'.b',X(1,label==2),X(2,label==2),'.r'); axis equal  %,X(1,label==3),X(2,label==3),'.g'); axis equal

% Generating the training date set
randn('seed',100)
n_of_points_per_group=100;

X_test=[mvnrnd(m11',S11,n_of_points_per_group); mvnrnd(m12',S12,n_of_points_per_group); mvnrnd(m21',S21,n_of_points_per_group)]'; %mvnrnd(m22',S22,n_of_points_per_group)]'; %mvnrnd(m3',S3,n_of_points_per_group)]';
label_test=[ones(1,n_of_points_per_group) ones(1,n_of_points_per_group) 2*ones(1,n_of_points_per_group)]; % 2*ones(1,n_of_points_per_group)];  % 3*ones(1,n_of_points_per_group)];
[l,p]=size(X_test);
%Plot the test data set
%figure; plot(X_test(1,label==1),X_test(2,label==1),'.b',X_test(1,label==2),X_test(2,label==2),'.r'); axis equal  %,X_test(1,label==3),X_test(2,label==3),'.g'); axis equal


y=ordinal(label);   %Converting the 'label' to ordinal values
y_test=ordinal(label_test);   %Converting the 'label' to ordinal values

no_of_base_classifiers=2000;

ens=fitensemble(X',y,'AdaBoostM1',no_of_base_classifiers,'Tree')

view(ens.Trained{8})

ens. ReasonForTermination

L=loss(ens,X',y,'mode','cumulative');

L_test=loss(ens,X_test',y_test,'mode','cumulative');

figure; plot(1:no_of_base_classifiers,L,'r',1:no_of_base_classifiers,L_test,'b')
 xlabel('Number of base classfiers');
     ylabel('Error');