%=====================================================================
% Exercise 8.39
% This script implements the CTA APSM algorithm
%=====================================================================

function [err, x]=APSM_distrib_sparse(inputvec,y,data,errfun)

 

L = data.L; N = data.N;
mu = data.mu;
if ~isa(mu,'function_handle')
    temp_mu = mu;
    mu = @(i) temp_mu;
end

nodes = data.nodes;

 
epsilon = data.epsilon;
h = data.h;

 q=data.q;

C = data.C; % Combination matrix

if ~isfield(data,'normalized')
    normalized = 0;
else
   normalized = data.normalized; 
end
if isfield(data,'initial_estimate'),
    for inod=1:nodes, x{inod} = data.initial_estimate; end
else for inod=1:nodes, x{inod} = zeros(L,1); end, end

a2v = zeros(1,nodes); 

for inod=1:nodes, err{inod} = zeros(1,N); end

%% transforms the input data to a function handle
if ~isa(inputvec,'function_handle')
    if length(inputvec)==N+L-1 % the input is shift invariant (system identification)
        inputseq = inputvec; inputvec = @(n) inputseq(L+n-1:-1:n); clear inputseq
    elseif size(inputvec,2)>1 % The imput matrix A is given
        inputmatrix = inputvec; inputvec = @(n) inputmatrix(:,n); clear inputmatrix
    end
end


for inod=1:nodes, y_n{inod}=[]; end

epsil = 0.1;
gradf = @(x) (1./(epsil+abs(x))).*sign(x);
f = @(x) sum(abs(x)./(epsil + abs(x)));
gradfmat = zeros(1,nodes);
fmat=zeros(1,nodes);

for i=1:N
     
       if mod(i,1000)==0, fprintf('%i, ',i),  end
     
   
      for inod = 1:nodes
        neighb = find(C(:,inod)~=0);
        
        x_c{inod} = zeros(L,1);
        for ii=1:length(neighb)%Combination Stage
            x_c{inod} = x_c{inod} + C(neighb(ii),inod)*x{neighb(ii)};
        end
    end
    for inod = 1:nodes
      
 
        
       
           Xq=inputvec(inod,max(1,i-q+1):i); 
           yvec=y{inod}(max(1,i-q+1):i);
           
           sx=length(max(1,i-q+1):i);
            sum1=zeros(L,1);
             sum2=0;
   
           for jj=1:sx
              
                p=projection(yvec(jj),Xq(:,jj), x_c{inod},epsilon);%Compute projection onto the hyperslab
     
             sum1=sum1+(1/sx)*p;
             sum2=sum2+(1/sx)*(((p-x_c{inod})'*(p-x_c{inod})));              
           end
          

    if sum1==x_c{inod}%Extrapolation Parameter computation
          Mn=1;
    else
      Mn=sum2/((sum1-x_c{inod})'*(sum1-x_c{inod}));          
    end

           
             
            x{inod}=x_c{inod}+0.2*Mn*(sum1-x_c{inod});%APSM update
          
            end

         
          
 

    %end
    
    
    if nargin(errfun)==2
        for inod=1:nodes
            err{inod}(i) = errfun(i,x{inod});
        end
    else
        for inod=1:nodes
            err{inod}(i) = errfun(x{inod});
        end
    end
    end
end

 

function [fnew]=projection(y,x,f,epsilon)%Function computing the projection of a vector onto a hyperslab
 
inner=f'*x;

if((inner-y<-epsilon))
  
    fnew=f+((y-inner-epsilon)/(x'*x))*x;
else if (inner-y>epsilon)
 %   err=f'*x-y;
 
    fnew=f+((y-inner+epsilon)/(x'*x))*x;
    else
  %      err=0;
        fnew=f;
    end
end
end