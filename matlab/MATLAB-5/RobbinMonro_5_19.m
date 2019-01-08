%===========================================================
% Exercise 5.19
% Plot the mean value of the obtaoned estimate  and the 
% corresponding variance 
% applying the Robins-Monro algorithm.
%==========================================================

close 
L=2;%Dimension of the unknown vector
N=500; %Number of Data
theta=randn(L,1);%Unknown parameter
w=zeros(L,1);%Initial Estimate 

IterNo=1000;
wtot=zeros(N,IterNo);
noisevar=0.1;
X = randn(L,N);%Input
inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);%Noise process

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;

for It=1:IterNo
    X = randn(L,N);
inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;
w=zeros(L,1);
        for i=1:N
          mu=1/(i); %Step size
        
            e=y(i)-w'*inputvec(i); %Error computation
            w=w+mu*e*inputvec(i);
        wtot(i,It)=w(1);
        end

end
theta1=theta(1)*ones(N,1);
plot(theta1,'r');hold on;
meanw=mean(wtot');
  plot(meanw,'k-');hold on
         for i=1:N
   
    if mod(i,10)==0 && i>30
    errorbar(i,meanw(i),std(wtot(i,:))) ;  hold on     
    end    
    end