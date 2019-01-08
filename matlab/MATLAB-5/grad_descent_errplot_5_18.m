%========================================================
% Exercise 5.18
% This function computes the error curves
% for the  Sigma=[1 0;0 0.1] covariance scenario.
% Two curves are presented, the first one
% corresponding to the optimum step-size \mu_o
% and the second one corresponding to \mu_o/2.
%==========================================================



function [xopt,fopt,niter,gnorm,dx] = grad_descent(varargin)
 close all
 clear
 
if nargin==0
    % define starting point
    x0 = [-2 -1]';
 elseif nargin==1
    % if a single input argument is provided, it is a user-defined starting
    % point.
    x0 = varargin{1};
else
    error('Incorrect number of input arguments.')
end

% termination tolerance
tol = 1e-8;

% maximum number of allowed iterations
maxiter = 200;

% minimum allowed perturbation
dxmin = 1e-8;


% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = x0; niter = 0; dx = inf;

 Sigma=[1 0;0 0.1]; % Covariance matrix. For the second part of the exercise set this covariance to Sigma=[1 0;0 1] 
 p=[0.05 .03]';

 
   alpha = 2/(1+0.1);%Compute the optimum stepsize

opt=double(inv(Sigma)*p);

f = @(x1,x2) (p(1)-Sigma(1,1)*x1)^2+(p(2)-Sigma(2,2)*x2)^2

colormap('gray');hold on

err=zeros(50,1);
 [X,Y] = meshgrid(-3:.2:3, -3:.2:3);
 
hold on
% gradient descent algorithm:
while  ( (niter <= maxiter ))
    % calculate gradient:
    err(niter+1)=double(norm(x-opt)^2);
    g = double(grad(x,Sigma,p));
     xnew = double(x - alpha*g);
     if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    refresh
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
     x = xnew;
     
    
end
xopt = x;
 niter = niter - 1;
  
plot(err,'b-')
hold on



Sigma=[1 0;0 0.1];
p=[0.05 .03]';

% define the objective function:
opt=double(inv(Sigma)*p);
% termination tolerance
tol = 1e-8;

% maximum number of allowed iterations
maxiter = 200;

% minimum allowed perturbation
dxmin = 1e-8;

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha =  1/(1+0.1);
   

% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = x0; niter = 0; dx = inf;

  opt=inv(Sigma)*p;
 
f = @(x1,x2) (p(1)-Sigma(1,1)*x1)^2+(p(2)-Sigma(2,2)*x2)^2
 
err=zeros(50,1);
 [X,Y] = meshgrid(-3:.2:3, -3:.2:3);
zz=(p(1)-Sigma(1,1)*X).^2+(p(2)-Sigma(2,2)*Y).^2;
sx=length(-3:.2:3);
 hold on
% gradient descent algorithm:
while   niter <= maxiter   
    % calculate gradient:
     err(niter+1)=double(norm(x-opt)^2);
    g = double(grad(x,Sigma,p));
    gnorm = norm(g);
      % take step:
    xnew = double(x - alpha*g);
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
     refresh
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
  dx=double(dx);
    x = double(xnew);
    xnew
   
    
end
xopt = x;
 niter = niter - 1;
  plot((err),'r-')

 % define the gradient of the objective
function g = grad(x,Sigma,p)
g=Sigma*x-p;
 