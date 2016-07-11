%=========================================================
%  Exercise 5.18
%  This function plots the convergence path together with the 
%  isovalue contours  of the  covariance matrix Sigma. 
%=========================================================
 
function [xopt,fopt,niter,gnorm,dx] = grad_descent(varargin)
close  all
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
maxiter = 1000;

% minimum allowed perturbation
dxmin = 1e-8;


% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = x0; niter = 0; dx = inf;


%Define the covariance matrix; 
Sigma=[1 0;0 0.1];
p=[0.05 .03]';

% Set the optimum stepsize;
   alpha =2/(0.1+1);

 f = @(x1,x2) (p(1)-Sigma(1,1)*x1)^2+(p(2)-Sigma(2,2)*x2)^2
 
colormap('gray');hold on


  [X,Y] = meshgrid(-3:.2:3, -3:.2:3);
zz=(p(1)-Sigma(1,1)*X).^2+(p(2)-Sigma(2,2)*Y).^2;
sx=length(-3:.2:3);
[c f]=contour(X,Y,zz,20);hold on
 while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = grad(x,Sigma,p);
    gnorm = norm(g);
     % take step:
    xnew = x - alpha*g;
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
  
    plot([x(1) xnew(1)],[x(2) xnew(2)],'r.-','LineWidth',1,'markers',15);hold on
     refresh
    % update termination metrics
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    xnew
end
niter
xopt = x;
 niter = niter - 1;

% define the gradient of the objective
function g = grad(x,Sigma,p)
g=Sigma*x-p;
 