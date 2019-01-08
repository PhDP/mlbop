%-----------------------------------------------------------------
%  Exercise 17.9
%  SIR particle filtering for various time instances.
%-----------------------------------------------------------------

clc; clear all; close all;

%% Initialization of Variables
N=200; % particle streams
x=zeros(1,N);
y=zeros(1,N);

%% State Model
x(1)=1;
for i=2:1:N
    x(i)=x(i-1)+randn;
    y(i)=randn+x(i);
end

%% Plot of State Model
figure(1)
hold all
plot(1:N,x, 'r')
plot(1:N,y, 'g')

%% Initialize Variables
par=200; %Number of Particles
ParPos=zeros(N,par); % Positions of Particles
weights=ParPos; % Weights of Particles
ParPos(1,:)=randn(1,par); % Generate random positions for particles
weights(1,:)=1/par; % weights of particles (initially uniformly distributed)
Xc=zeros(1,par); %Vector to store Resampling Positions

Predict=zeros(1,N); %Prediction Vector
Predict(1)=ParPos(1,:)*weights(1,:)';

%% SIR algorithm
for i=2:N
    for j=1:par
        ParPos(i,j)=normrnd(x(i)-ParPos(i-1,j),1); %Position Update
        weights(i,j)=weights(i-1,j)*exp(-((y(i)-ParPos(i,j))^2)/2); %weight update
    end
    weights(i,:)=weights(i,:)/sum(weights(i,:)); %Normalization of Weights
    Meff=1/(weights(i,:)*weights(i,:)'); %effective Particle Size
    Predict(i)=ParPos(i,:)*weights(i,:)'; %Prediction
    %% Resampling
    if Meff<par/2
        for d=1:par
            km=mnrnd(1,weights(i,:),1); %Drawing from a multivariable distribution
            Xc(d)=ParPos(i,:)*km'; %Select the Position that is drawn
        end
        %% Update Positions and Weights
        ParPos(i,:)=Xc;
        weights(i,:)=1/par;
    end
end

for fig=[2:4, 30, N]
    figure(fig)
    stem(1:par, weights(fig, :))
end


figure(1)
plot(1:N,Predict,'--k')