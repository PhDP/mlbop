%=====================================================================
% Exercise 8.38
% This script generates the model, the data and runs 
% the APSM, the RLS, the APA and the NLMS algorithms.  
%=====================================================================

function  APSM 
 figure 
L=200;%Dimension of the unknown vector
N=3500; %Number of Data
theta=randn(L,1);%Unknown parameter
 

IterNo=100;
 MSE1=zeros(N,IterNo);
MSE2=zeros(N,IterNo);
MSE3=zeros(N,IterNo);
MSE4=zeros(N,IterNo);


noisevar=0.01;%For the high noise scenario set this value equal to 0.3
epsilon=sqrt(2)*sqrt(noisevar); %Epsilon parameter used in the APSM

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

X = convmtx(xcorrel,L)';%Generate noise process
X(:,1:L-1) = [];

inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);%Generate noise process

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;
w=zeros(L,1);%Estimate vector
   mu=0.2;
   delta=0.001;
   q=30;% Number of window used for the APA and the APSM.  
        for i=1:N
         
            if i>q
            qq=i:-1:i-q+1;
            
            yvec=y(qq);
            Xq=inputvec(qq);
           Xq=Xq'; 
          
            e=yvec-Xq*w;
           eins=y(i)-w'*inputvec(i);
            w=w+mu*Xq'*inv(delta*eye(q)+Xq*Xq')*e;%APA update
            MSE1(i,It)=eins^2;
            end
        end
        w=zeros(L,1);
  
   for i=1:N
         
            if i>q
            qq=i:-1:i-q+1;
            yvec=y(qq);
            Xq=inputvec(qq);
            Xq=Xq';  
           
         sum1=zeros(L,1);
         sum2=0;

           for jj=1:q
              
                p=projection(yvec(jj),Xq(jj,:),w,epsilon);%Compute projection onto hyperslab
     
             sum1=sum1+(1/q)*p;
             sum2=sum2+(1/q)*(((p-w)'*(p-w)));              
           end
          

    if sum1==w%Extrapolation Parameter Computation
          Mn=1;
    else
      Mn=sum2/((sum1-w)'*(sum1-w));          
        end

           
            eins=y(i)-w'*inputvec(i);
            w=w+Mn*0.5*(sum1-w);%APSM update
            MSE2(i,It)=eins^2;
            end
   end
        
        
    w=zeros(L,1);%RLS recursion
        delta=0.001;
   Pi=delta*eye(L);
   P=(1/delta)*eye(L);
        for i=1:N

            gamma=1/(1+inputvec(i)'*P*inputvec(i));
            gi=P*inputvec(i)*gamma;
           e=y(i)-w'*inputvec(i);
            w=w+gi*(e);%RLS update
            P=P-(gi*gi')/gamma;
            MSE3(i,It)=e^2;
        end
        
        w=zeros(L,1);%NLMS Recursion
        delta=0.001;
         mu=1.2;
        for i=1:N
          
            e=y(i)-w'*inputvec(i);
           mun=mu/(delta+inputvec(i)'*inputvec(i));
            w=w+mun*e*inputvec(i);%NLMS update
            MSE4(i,It)=e^2;
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

function [fnew]=projection(y,x,f,epsilon)%Projection onto the hyperslab function

 x=x';
inner=f'*x;
if((inner-y<-epsilon))

    fnew=f+((y-inner-epsilon)/(x'*x))*x;
else if (inner-y>epsilon)
  
    fnew=f+((y-inner+epsilon)/(x'*x))*x;
    else
         fnew=f;
    end
end

