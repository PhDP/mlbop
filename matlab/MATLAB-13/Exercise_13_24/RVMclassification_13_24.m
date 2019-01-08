%-----------------------------------------------------------------
%  Exercise 13.24
%  RVM classification
%  Use SB1_Release_110 Matlab package
%      found at http://www.miketipping.com/sparsebayes.htm.
%-----------------------------------------------------------------


clc; clear; close all; format compact; format long eng;

rng('default');

% data dimension
l = 2; 

% number of data points
N = 150; 

% generate the data
x1 = 10 * rand(l,N) - 5; 
y1 = zeros(N,1);
for i = 1 : N 
    t = .05 * (x1(1,i)^3 + x1(1,i)^2 + x1(1,i) + 1);
    if t + 2*randn(1) > x1(2,i)
        y1(i) = 1 ;
    else
        y1(i) = 0 ;
    end
end

% Plot the training data 
COL_data1	= 0.5*ones(1,3);
COL_data2	= [1 0 0];
COL_boundary25	= 0.5*ones(1,3);
COL_boundary75	= 0.5*ones(1,3);
COL_rv		= 'k';
figure(1)
whitebg(1,'w')
clf
h_c1	= plot(x1(1,y1==0),x1(2,y1==0),'.','MarkerSize',10,'Color',COL_data1);
hold on
h_c2	= plot(x1(1,y1==1),x1(2,y1==1),'.','MarkerSize',10,'Color',COL_data2);
box	= 1.1*[min(x1(1,:)) max(x1(1,:)) min(x1(2,:)) max(x1(2,:))];axis(box)
set(gca,'FontSize',12)
drawnow


% Set verbosity of output (0 to 4)
setEnvironment('Diagnostic','verbosity',3);
% Set file ID to write to (1 = stdout)
setEnvironment('Diagnostic','fid',1);

kernel_	= 'gauss';
% variance of the gaussian kernel
width		= 3; 
% max iterations
maxIts	= 1000; 

% Set up initial hyperparameters - precise settings should not be critical
initAlpha	= (1/N)^2;
% Set beta to zero for classification
initBeta	= 0;
useBias	= true;
monIts		= round(maxIts/10);

% rvm classifier
[weights, used, bias, marginal, alpha, beta, gamma] = SB1_RVM(x1.',y1,initAlpha,initBeta,kernel_,width,useBias,maxIts,monIts);

% visualise the results over a grid
gsteps		= 50;
range1		= box(1):(box(2)-box(1))/(gsteps-1):box(2);
range2		= box(3):(box(4)-box(3))/(gsteps-1):box(4);
[grid1, grid2]	= meshgrid(range1,range2);
Xgrid		= [grid1(:) grid2(:)];

% Evaluate RVM
PHI		= SB1_KernelFunction(Xgrid,x1(:,used).',kernel_,width);
y_grid		= PHI*weights + bias;

% apply sigmoid for probabilities
p_grid		= 1./(1+exp(-y_grid)); 

% show decision boundary (p=0.5) and illustrate p=0.25 and 0.75
[~,h05]		= contour(range1,range2,reshape(p_grid,size(grid1)),[0.5],'-');
set(h05, 'Color',COL_rv,'LineWidth',1);

% show relevance vectors
h_rv	= plot(x1(1,used),x1(2,used),'o','LineWidth',2,'MarkerSize',10, 'Color','k');

hold off

