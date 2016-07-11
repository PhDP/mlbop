%===================================================================
% Exercise 8.40
% Implementation of the Pegasos Algorithm.
%==================================================================

function [wT,b] = pegasos(X,Y,lamda,k,maxIter,Tolerance)
% X imput matrix
% Y labels
% lamda stepsize parameter of pegasos
% k data processed per iteration step
% maxIter the maximum number of iterations 
% Tolerance Error Tolerance 
[N,d]=size(X);
size(Y,1)
 
if(size(Y,1)~=N)
    fprintf('\nError: Number of elements in X and Y must same\nSee pegasos usage for further help\n');
    return;
end
if(sum(Y~=1 & Y~=-1)>0)
    fprintf('\nError: Y must be 1 or -1\nSee pegasos usage for further help\n');
    return;
end

if(nargin<3 || isempty(lamda)),    lamda=1;  end
if(nargin<4 || isempty(maxIter)),    maxIter=10000;  end
if(nargin<5 || isempty(k)),    k=ceil(0.1*N);  end
if(nargin<6 || isempty(Tolerance)),    Tolerance=10^-6;  end

w=rand(1,size(X,2));
w=w/(sqrt(lamda)*norm(w));%Initialization
for t=1:maxIter
     b=mean(Y-X*w(t,:)');
    idx=randint(k,1,[1,size(X,1)]);
    At=X(idx,:);
    yt=Y(idx,:);
    idx1=(At*w(t,:)'+b).*yt<1;
    etat=1/(lamda*t);
    w1=(1-etat*lamda)*w(t,:)+(etat/k)*sum(At(idx1,:).*repmat(yt(idx1,:),1,size(At,2)),1);
    w(t+1,:)=min(1,1/(sqrt(lamda)*norm(w1)))*w1;%Updated estimate for the PEGASOS algorithm
    if(norm(w(t+1,:)-w(t,:)) < Tolerance)
        break;
    end
end
if(t<maxIter)
    fprintf('\nW converged in %d iterations.',t);
else
    fprintf('\nW not converged in %d iterations.',maxIter);
end
wT=mean(w,1);
% wT=w(end,:);
b=mean(Y-X*wT');
Tr=sum(sign(X*wT'+b)==Y);
F=size(X,1)-Tr;
TrainAccuracy=100*Tr/(Tr+F);
fprintf('\nPegasos Accuracy on Training set = %.4f %%\n',TrainAccuracy);
end
