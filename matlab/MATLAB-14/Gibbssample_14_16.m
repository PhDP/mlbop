%-----------------------------------------------------------------
%  Exercise 14.16
%  Gibbs sampling
%-----------------------------------------------------------------

rand('seed' ,12345);
n_samples = 10^3;
 
MeanValue = [0 0]; % mean vector
sigma(1) = 0.5; 
sigma(2) = 0.5; 
 
propSigma = 1; % proposal variance
minimum = [-4 -4];
maximum = [4 4];
 
%% Samples initialization
x = zeros(n_samples,2);
x(1,1) = unifrnd(minimum(1), maximum(1));
x(1,2) = unifrnd(minimum(2), maximum(2));
 
dimensions = 1:2; 

%% Start Gibbs Sampler
t = 1;
while t < n_samples
    t = t + 1;
    T = [t-1,t];
    for iD = 1:2 % 
        % sample update
        nIx = dimensions~=iD;  
        % conditional mean
        ConditionalMean = MeanValue(iD) + sigma(iD)*(x(T(iD),nIx)-MeanValue(nIx));
        % conditional variance
        ConditionalVariance = sqrt(1-sigma(iD)^2);
        % sample from conditional pdf
        x(t,iD) = normrnd(ConditionalMean,ConditionalVariance);
    end
end

%% Plot Samples
figure(1);
h1 = scatter(x(:,1),x(:,2),'r.');

%% Plot Figures
hold on;
for t = 1:10 %change limit to plot the number of samples in h2 legend
    plot([x(t,1),x(t+1,1)],[x(t,2),x(t,2)],'k-');
    plot([x(t+1,1),x(t+1,1)],[x(t,2),x(t+1,2)],'k-');
    h2 = plot(x(t+1,1),x(t+1,2),'ko');
end
 
legend([h1,h2],{'Samples','First 10 Samples'},'Location','Northwest')
xlabel('x_1');
ylabel('x_2');
grid on