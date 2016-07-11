%-----------------------------------------------------------------
%  Exercise 17.8
%  SIS particle filtering algorithm
%-----------------------------------------------------------------

clc; clear all; close all;

%% Initialization of Variables
N=90; % stream of particles
x=zeros(1,N);
y=x;

%% State Model
x(1)=1;
y(1)=x(1)+randn;

for i=2:N
    x(i)=x(i-1)+randn(1);
    y(i)=x(i)+randn(1);
    
end

%% Figures of State Model
figure(1)
plot(1:N,x)
figure(2)
plot(1:N,y)

%% Main Algorithm
par=1000;
ParPos=randn(1,par);
weights(1:par)=1/par;

% figure(1) % Ucomment to plot the initial weights
% stem(ParPos,weights)

%% Plot of weights for different time snapshots
for i=[2:4, 10, 30] % Change i to see the weights for the
    % desired time instances.
    for j=1:par
        ParPos(j)=abs(x(i)-ParPos(j))+randn;
        weights(j)=weights(j)*exp(-((y(i)-ParPos(j))^2)/2);
    end
    weights=weights/sum(weights);
    figure(i) % generate a figure for each shanpshot of the algorithm
    stem(ParPos,weights)
end

