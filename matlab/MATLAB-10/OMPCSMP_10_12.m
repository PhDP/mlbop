%-----------------------------------------------------------------
% Exercise 10.12
% Implements the OMP and CSMP algorithms.
%----------------------------------------------------------------


function Matlab10_12

% Chooze a fixed signal length.
l = 100;


alpha = 0.2; % set alpha equal to 0.8 for question (b)
N = round(alpha*l);
error = zeros(2,10);
beta = 0.1:0.1:1;

X = randn(N,l)*sqrt(1/N);
for count = 1:10
   k = round(beta(count)*N);
   t = 2*k;

   % Make the unknown sparse vector
   theta = zeros(l,1);
   theta(randperm(l,k))=randn(k,1);

       
   y = X*theta;
       
%  Uncomment the next 4 lines in order to add some noise (question (c))
%        SNR = 20;
%        varnoise = (k/N)/10^(SNR/10);
%        noise = randn(size(y))*sqrt(varnoise);
%        y = y + noise;
   
       % Run OMP
       theta_OMP = OMP(X,y,k);
       error(1,count) = norm(theta - theta_OMP);
       
       % Run CSMP
       theta_CSMP = CSMP(X,y,k,t);
       error(2,count) = norm(theta - theta_CSMP);

end
   
figure;plot(beta,error(1,:),beta,error(2,:))
legend('OMP','CSMP')
end




function theta=CSMP(X,y,k,t)
   
l = size(X,2);
if t>l, t=l; end % check if t is larger than l


    theta = zeros(l,1);
    residual=y;
  %  normx = sqrt(sum(X.^2,1))';
    for i=1:1:20,
        b1 = find(theta);
        proxy = X'*residual;
        [~,b2] = sort(abs(proxy),'descend');
        S = unique([b1;b2(1:t)]);
        
        theta_=pinv(X(:,S))*y;

        theta = zeros(size(X,2),1);
        theta(S) = theta_;
        [~,b] = sort(abs(theta),'descend');
        theta(b(k+1:end)) = 0;
        residual=y-X*theta;
        
    end;
end

function theta=OMP(X,y,k)

    residual=y;
    S=zeros(k,1);
    normx = sqrt(sum(X.^2,1))';
    for i=1:1:k,
        proj=X'*residual;
        proj = proj./normx;
        [~,pos]=max(abs(proj));
        pos=pos(1);
        S(i)=pos;
        theta_=pinv(X(:,S(1:i)))*y;

        theta = zeros(size(X,2),1);
        theta(S(1:i)) = theta_;
        residual=y-X*theta;
        
      %  residual=y-X(:,indx(1:i))*theta_; %a faster implementation
    end;
   % theta = zeros(size(X,2),1);
% theta(indx) = theta_;
end
