%=====================================================================
% Exercise 5.22
% This function generates the model  and runs the LMS, the NLMS
% and the APA algorithm for two different choices of the parameter q.
%======================================================================


function   APA
  figure 
L=60;%Dimension of the unknown vector
N=3500; %Number of Data
theta=randn(L,1);%Unknown parameter
w=zeros(L,1);%Estimate Vector

IterNo=100;
 MSE1=zeros(N,IterNo);
MSE2=zeros(N,IterNo);
MSE3=zeros(N,IterNo);
MSE4=zeros(N,IterNo);


noisevar=0.01;
X = randn(L,N);
inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;

for It=1:IterNo
    xc = randn(N+L-1,1);
    xc = xc/std(xc);

X = convmtx(xc,L)';
X(:,1:L-1) = [];

inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;
w=zeros(L,1);
   mu=0.1;
   delta=0.001;
   q=30;
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
   q=10;
        for i=1:N
         
            if i>q
            qq=i:-1:i-q+1;
            
            yvec=y(qq);
            Xq=inputvec(qq);
           Xq=Xq'; 
          
            e=yvec-Xq*w;
          eins=y(i)-w'*inputvec(i);
            w=w+mu*Xq'*inv(delta*eye(q)+Xq*Xq')*e;
            MSE2(i,It)=eins^2;
            end
        end
   
        
        
        w=zeros(L,1);
        delta=0.001;
         mu=0.35;
        for i=1:N
         
          %  mu=1;%/(i^0.5);
            e=y(i)-w'*inputvec(i);
           mun=mu/(delta+inputvec(i)'*inputvec(i));
            w=w+mun*e*inputvec(i);
            MSE4(i,It)=e^2;
        end

          w=zeros(L,1);
         mu=0.025;
        for i=1:N
         
          %  mu=1;%/(i^0.5);
            e=y(i)-w'*inputvec(i);
            w=w+mu*e*inputvec(i);
            MSE3(i,It)=e^2;
        end

        
end
MSEav1=sum(MSE1')/IterNo;
MSEav2=sum(MSE2')/IterNo;
MSEav3=sum(MSE3')/IterNo;
MSEav4=sum(MSE4')/IterNo;

plot(10*log10(MSEav1),'r');hold on
plot(10*log10(MSEav2),'g');hold on
plot(10*log10(MSEav3),'b');hold on
plot(10*log10(MSEav4),'k')

