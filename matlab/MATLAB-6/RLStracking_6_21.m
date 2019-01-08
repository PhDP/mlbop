%====================================================================
% Exercise 6.21
% Time-varying parameter estimation using the RLS and the
% NLMS.
%====================================================================

 
function RLStrack
 close all
 clear
L=5;%Dimension of the unknown vector
N=1000; %Number of Data
theta=zeros(L,N);%Unknown parameter
  

IterNo=200;

MSE1=zeros(N,IterNo);
MSE2=zeros(N,IterNo);


ara=0.97 ;
qqm=0.1; %Parameters for the time varying channel



noisevar=0.001;
epsilon=sqrt(2)*noisevar;



for It=1:IterNo
    It
    theta=zeros(L,N);
    theta(:,1)=randn(L,1);
for ii=2:N%Generate the time varying channel
    theta(:,ii)=ara*theta(:,ii-1)+qqm*randn(L,1);
end
    xcorrel = randn(N+L-1,1);
 
xcorrel = xcorrel/std(xcorrel);
% 
 
X=randn(L,N);
inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
for ii=1:N

y(ii) = X(:,ii)'*theta(:,ii);

end

y = y + noise;
     w=zeros(L,1);
        delta=0.001;
        lambda=0.995;
  
   P=(1/delta)*eye(L);
        for i=1:N
 
          P=(1/lambda)*(P- ((1/lambda)*P*inputvec(i)*inputvec(i)'*P)/(1+(1/lambda)*inputvec(i)'*P*inputvec(i)));  
          e=y(i)-w'*inputvec(i);
            w=w+P*inputvec(i)*(e);
     
            MSE1(i,It)=e^2;
        end
        
        w=zeros(L,1);
        delta=0.001;
         mu=0.5;
        for i=1:N
         
             e=y(i)-w'*inputvec(i);
           mun=mu/(delta+inputvec(i)'*inputvec(i));
            w=w+mun*e*inputvec(i);
            MSE2(i,It)=e^2;
        end
 

        
end
  MSEav1=sum(MSE1')/IterNo;
 MSEav2=sum(MSE2')/IterNo;
 plot(10*log10(MSEav1),'b');hold on
plot(10*log10(MSEav2),'k')
 