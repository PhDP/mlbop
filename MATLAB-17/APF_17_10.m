%-----------------------------------------------------------------
%  Exercise 17.10
%  Auxiiliary particle filtering.
%-----------------------------------------------------------------

clc; clear all; close all;
%% Initialization of Variables
N=50; % time instances
x=zeros(1,N);
y=x;
x(1)=1;
y(1)=randn+x(1);

%% State Model
for i=2:N
    x(i)=x(i-1)+randn;
    y(i)=x(i)+randn;
end

%% Figures of State Model 
% figure(1) % Uncomment to plot the state model
% hold on
% plot(1:N,x,'g')
% plot(1:N,y,'r')
% hold off

%% Initialize Variables
pd='Normal';
par=2000; %Number of Particles
ParPos=zeros(N,par); % Positions of Particles
weights=ParPos; % Weights of Particles
ParPos(1,:)=randn(1,par)+x(1);
weights(1,:)=1/par;
xHat=zeros(1,par);
kweights=zeros(1,par);
xHat2=xHat;

Prediction=zeros(1,N);
Prediction(1)=ParPos(1,:)*weights(1,:)';

%% Main APF algorithm
for i=2:N
    xHat=x(i)-ParPos(i-1,:); % Expectation of x_n given previous particles
    kweights(1:par)=weights(i-1,:).*pdf(pd,y(i)-xHat,1);
    kweights=kweights/sum(kweights);
    for j=1:par
        Km=mnrnd(1,kweights,1);
        ParPos(i,j)=randn+x(i)-ParPos(i-1,:)*Km';
        xHat2(j)=xHat*Km';
    end
    weights(i,:)=pdf(pd,y(i)-ParPos(i,:),1)./pdf(pd,y(i)-xHat2,1);
    weights(i,:)=weights(i,:)/sum(weights(i,:));
    Prediction(i)=ParPos(i,:)*weights(i,:)';
end
figure(2)
stem(ParPos(2,:),weights(2,:),'r')
figure(4)
stem(ParPos(4,:),weights(4,:),'r')
figure(10)
stem(ParPos(10,:),weights(10,:),'r')
figure(30)
stem(ParPos(30,:),weights(30,:),'r')

