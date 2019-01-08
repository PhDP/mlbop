%=======================================================================
% Exercise 6.20
% This function generates the model, the data and runs the RLS, the NLMS
% and the APA algorithms.
%========================================================================

function RLS
figure
L=200;%Dimension of the unknown vector
N=3500; %Number of Data
theta=randn(L,1);%Unknown parameter
 

IterNo=100;
 MSE1=zeros(N,IterNo);
MSE2=zeros(N,IterNo);
MSE3=zeros(N,IterNo);
 

noisevar=0.01;
epsilon=sqrt(2)*noisevar;

X = randn(L,N);
inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;

for It=1:IterNo
    It
    xcorrel = randn(N+L-1,1);
xcorrel = xcorrel/std(xcorrel);

X = convmtx(xcorrel,L)';
X(:,1:L-1) = [];

inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;
w=zeros(L,1);
   mu=0.2;
   delta=0.001;
   q=30;% Number of window used for the APA
        for i=1:N
         
            if i>q
            qq=i:-1:i-q+1;
            
            yvec=y(qq);
            Xq=inputvec(qq);
           Xq=Xq'; 
      
         
            e=yvec-Xq*w;
           eins=y(i)-w'*inputvec(i);
            w=w+mu*Xq'*inv(delta*eye(q)+Xq*Xq')*e;
            MSE1(i,It)=eins^2;
            end
        end
        w=zeros(L,1);
  
           
        
    w=zeros(L,1);%RLS recursion
        delta=0.001;
    P=(1/delta)*eye(L);
        for i=1:N

            gamma=1/(1+inputvec(i)'*P*inputvec(i));
            gi=P*inputvec(i)*gamma;
            e=y(i)-w'*inputvec(i);
            w=w+gi*(e);
            P=P-(gi*gi')/gamma;
            MSE2(i,It)=e^2;
        end
        
        w=zeros(L,1);%NLMS Recursion
        delta=0.001;
         mu=1.2;
        for i=1:N
         
          %  mu=1;%/(i^0.5);
            e=y(i)-w'*inputvec(i);
           mun=mu/(delta+inputvec(i)'*inputvec(i));
            w=w+mun*e*inputvec(i);
            MSE3(i,It)=e^2;
        end
 

        
end
MSEav1=sum(MSE1')/IterNo;
MSEav2=sum(MSE2')/IterNo;
MSEav3=sum(MSE3')/IterNo;
 
plot(10*log10(MSEav1),'r');hold on
plot(10*log10(MSEav2),'g');hold on
plot(10*log10(MSEav3),'b');hold on
 
 