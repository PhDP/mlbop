%-----------------------------------------------------------------
%  Exercise 17.11
%  Stochastic volatility model.
%-----------------------------------------------------------------

clc; clear all; close all;
%% Initialization of Variables
N=50; % time instances
x=zeros(1,N);
y=zeros(1,N);



%% Stochastic Volatility Model
a = 0.97;
b = 0.69;
sigma_h = sqrt(0.178);
sigma_v =1; 
x(1)=1;
y(1)=b*(sigma_v*randn)*exp(x(1)/2);
for i=2:N
    x(i)=a*x(i-1)+(sigma_h*randn);
    y(i)=b*(sigma_v*randn)*exp(x(i)/2);
end

%% Plot Figures of Input and Output Sequences
figure(1)
hold on
%plot(1:N,x,'g')
plot(1:N,y,'r')
hold off

%% APF Algorithm
pd='Normal';
par=1000; %Number of Particles
ParPos=zeros(N,par); % Positions of Particles
weights=ParPos; % Weights of Particles
ParPos(1,:)=randn(1,par)+x(1);
weights(1,:)=1/par;
xHat=zeros(1,par);
kweights=zeros(1,par);
xHat2=xHat;

Prediction=zeros(1,N);
Prediction(1)=ParPos(1,:)*weights(1,:)';

for i=2:N
    xHat=x(i)-ParPos(i-1,:); %Expectation of x_n given previous particles
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

figure(1)
hold on
plot(1:N,Prediction,'--k')
grid on
%legend('y','x','\hat{x}')
hold off