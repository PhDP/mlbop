function [err, x]=LMS_distrib_sparse(inputvec,y,data,errfun)

 
L = data.L; N = data.N;
mu = data.mu;
if ~isa(mu,'function_handle')
    temp_mu = mu;
    mu = @(i) temp_mu;
end

nodes = data.nodes;

gamma = data.gamma;
 
h = data.h;

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
     
    if gamma==0
        for inod = 1:nodes
            gradfmat(inod) = norm(gradf(x{inod}))^2;
            fmat(inod) = f(x{inod});
        end
    end
    
    for inod = 1:nodes
      
         a_n=inputvec(inod,i); 
      
           er = (y{inod}(i) - a_n'*x{inod});

if normalized%adaptation stage. IF normalized is true the stepsize is normalized.
    x_c{inod} = x{inod} + mu(i)/norm(a_n)^2 * a_n * er ;
else
    x_c{inod} = x{inod} + mu(i) * a_n * er ;
end

    end
    for inod = 1:nodes %Combination stage
        neighb = find(C(:,inod)~=0);
        
        x{inod} = zeros(L,1);
        for ii=1:length(neighb)
            x{inod} = x{inod} + C(neighb(ii),inod)*x_c{neighb(ii)};
        end
    end
    
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


%fprintf('\n')
end
