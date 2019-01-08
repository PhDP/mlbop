%================================================================
% Exercise 5.21.
% This function generates the model 
% and runs the LMS and the TDLMS algorithms 
% for the AR process input scenario. 
%=================================================================

function TDLMS
  close 
L=10;%Dimension of the unknown vector
N=1500; %Number of Data
theta=randn(L,1);%Unknown parameter
w=zeros(L,1);%Estimate vector

IterNo=100;
 
MSE1=zeros(N,IterNo);
MSE2=zeros(N,IterNo);
noisevar=0.01;
correlcoeff=0.85;

delta=0.01;
beta=0.5;
T=dctmtx(L);%DCT matrix to be used for the TDLMS
for It=1:IterNo
  xcorrel = randn(N+L-1,1);
xcorrel = filter(1,[1 correlcoeff],xcorrel); % AR
xcorrel = xcorrel/std(xcorrel);

X = convmtx(xcorrel,L)';
X(:,1:L-1) = [];

inputvec = @(n) X(:,n);


inputvec = @(n) X(:,n);

noise = randn(N,1)*sqrt(noisevar);

y = zeros(N,1);
y(1:N) = X(:,1:N)'*theta;
y = y + noise;
w=zeros(L,1);
   mu=0.01;
   
        for i=1:N
         
          
            e=y(i)-w'*inputvec(i);
            w=w+mu*e*inputvec(i);
            MSE1(i,It)=e^2;
        end
        
        w=zeros(L,1);
         mu=0.01;
         sigmasquare=delta*ones(L,1);
   Dinv=eye(L);
         for i=1:N
            intran=T'*inputvec(i);
            e=y(i)-w'*intran;
            w=w+mu*e*Dinv*intran;%Compute the updated estimate
            sigmasquare = beta* sigmasquare+(1-beta)*abs(intran).^2;%Update the weights used in the TDLMS
            Dinv=diag(1./sigmasquare);
            MSE2(i,It)=e^2;
        end

end
MSEav1=sum(MSE1')/IterNo;
MSEav2=sum(MSE2')/IterNo;

plot(10*log10(MSEav1),'r');hold on
plot(10*log10(MSEav2),'g')
