%-----------------------------------------------------------------
%  Exercise 13.26
%  nonparametric Gaussian mixture example
%  Use the vdpgm package available at 
%   http://sites.google.com/site/kenichikurihara/academic-software
%-----------------------------------------------------------------


clear; clc; format long eng; format compact;

rng('default');

% 2-d data
D = 2; 
Ni = randi(30,5,1);

% # data points
N = sum(Ni); 

% # of gaussians
K = 5; 

% generate the data from a DP Gaussian mixture
mu = [-12.5 -4 2 10 3;2.5 -.1 -3.5 8 3];
Sigma1 = [1.4 .81;.81 1.3]; Sigma2 = [1.5 .2;.2 2.1]; Sigma3 = [1.6 1;1 2.9]; Sigma4 = [.5 .22;.22 .8]; Sigma5 = [1.5 1.4;1.4 2.4];
x1 = mvnrnd(mu(:,1), Sigma1, Ni(1));
x2 = mvnrnd(mu(:,2), Sigma2, Ni(2));
x3 = mvnrnd(mu(:,3), Sigma3, Ni(3));
x4 = mvnrnd(mu(:,4), Sigma4, Ni(4));
x5 = mvnrnd(mu(:,5), Sigma5, Ni(5));
% plot the data
% figure;plot(x1(:,1),x1(:,2),'ko',x2(:,1),x2(:,2),'ko',x3(:,1),x3(:,2),'ko',x4(:,1),x4(:,2),'ko',x5(:,1),x5(:,2),'ko');

X = [x1.' x2.' x3.' x4.' x5.' ];
% plot(X.','o')

% set parameters
[test_data1, test_data2] = meshgrid(-20:.1:15, -8:.1:12);
test_data = [test_data1(:).'; test_data2(:).'];
opts = mkopts_bj(10);
% set number of iterations
opts.sis = 1;
opts.get_q_of_z = 1;
% prediction results using 1 iteration
results = vdpgm(X, opts);
results.test_data = test_data;
results_predictive = vdpgm(X, results);
ppdf = reshape(results_predictive.predictive_posterior(:),size(test_data1,1),size(test_data1,2));

% plot the results
figure;
axh = gca;
set(axh,'XGrid','on','XColor', [1 1 1]);
set(axh,'YGrid','on');
C = repmat(flipdim(min(1,abs(1.2 - logspace(0.3,1,256)/10)).',1), 1,3);
colormap(C);
contourf(-20:.1:15, -8:.1:12, ppdf, 'Linestyle','none'); grid on;
hold on; 
plot(x1(:,1),x1(:,2),'ro',x2(:,1),x2(:,2),'ro',x3(:,1),x3(:,2),'ro',x4(:,1),x4(:,2),'ro',x5(:,1),x5(:,2),'ro');
hold off;


% set number of iterations
opts.sis = 2;
% prediction results using 2 iterations 
results = vdpgm(X, opts);
results.test_data = test_data;
results_predictive = vdpgm(X, results);
ppdf = reshape(results_predictive.predictive_posterior(:),size(test_data1,1),size(test_data1,2));

% plot the results
figure;
axh = gca;
set(axh,'XGrid','on','XColor', [1 1 1]);
set(axh,'YGrid','on');
C = repmat(flipdim(min(1,abs(1.2 - logspace(0.3,1,256)/10)).',1), 1,3);
colormap(C);
contourf(-20:.1:15, -8:.1:12, ppdf, 'Linestyle','none'); grid on;
hold on; 
plot(x1(:,1),x1(:,2),'ro',x2(:,1),x2(:,2),'ro',x3(:,1),x3(:,2),'ro',x4(:,1),x4(:,2),'ro',x5(:,1),x5(:,2),'ro');
hold off;

% set number of iterations
opts.sis = 5;
% prediction results using 5 iterations 
results = vdpgm(X, opts);
results.test_data = test_data;
results_predictive = vdpgm(X, results);
ppdf = reshape(results_predictive.predictive_posterior(:),size(test_data1,1),size(test_data1,2));

% plot the results
figure;
axh = gca;
set(axh,'XGrid','on','XColor', [1 1 1]);
set(axh,'YGrid','on');
C = repmat(flipdim(min(1,abs(1.2 - logspace(0.3,1,256)/10)).',1), 1,3);
colormap(C);
contourf(-20:.1:15, -8:.1:12, ppdf, 'Linestyle','none'); grid on;
hold on; 
plot(x1(:,1),x1(:,2),'ro',x2(:,1),x2(:,2),'ro',x3(:,1),x3(:,2),'ro',x4(:,1),x4(:,2),'ro',x5(:,1),x5(:,2),'ro');
hold off;
