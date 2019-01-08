
%---------------------------------------------------------------------
% Eercise 11-21
% Support Vector Machine classification
%
% This file uses a simple SMO implementation to train the SVM.
%---------------------------------------------------------------------


clear all;

rng(0);


disp('Starting');
N=150;

tic;

%---------Generate Samples-------------------------
f = @(x) (0.05*x.^3 + 0.05*x.^2 + 0.05*x + 0.05);

X_ = 10*rand(2,N) - 5;
x_ = X_(1,:);
y_ = X_(2,:);
fx = f(x_);



C1=[];
C2=[];
y_training = zeros(N,1);
for i=1:N 
    if y_(i) < fx(i) + 2*randn(1)
        C1 = [C1; i];
        y_training(i) = +1;
    else
        C2 = [C2; i];
        y_training(i) = -1;
    end;
end;

x_training = [x_;
              y_]';

%-----------------------------------------------------
          
          
%FOR N = 150
%Experiment1 C = 20
%Experiment2 C = 1

%FOR N=500
%Experiment1  C=1
%Experiment1  C=20


%-----Training---------
C=1;   
epsilon=0.01;
kernel_type='gaus';
kernel_params=sqrt(100);

[a, b] = SMO_classification(x_training, y_training, C, epsilon, kernel_type, kernel_params);
%-----------------------

%classifier
PP=1.1;
leftX = PP*min(x_);
rightX = PP*max(x_);
leftY = PP*min(y_);
rightY = PP*max(y_);
step=0.5;
[X,Y] = meshgrid(leftX:step:rightX, leftY:step:rightY);
[K,L]=size(X);
Z=zeros(K,L);
for i=1:K
    for j=1:L
        sum=0;
        for n=1:N
            sum = sum + a(n)*y_training(n)*kappa(x_training(n,:), [X(i,j) Y(i,j)], kernel_type, kernel_params);
        end;
        Z(i,j) = sum-b;
    end;
end;

%norm of solution
norm=0;
for n=1:N
    for m=1:N
        norm = norm + a(n)*a(m)*y_training(n)*y_training(m)*kappa(x_training(n,:), x_training(m,:), kernel_type, kernel_params);
    end;
end;






%find support vectors
SV = [];
for n=1:N
    if abs(a(n))>0.001
        SV = [SV n];
    end;
end;

SV1 = intersect(C1, SV);
SV2 = intersect(C2, SV);




figure(1);
hold on;

%plot classifier
v=[-1 0 1];
[Co,h] = contour(X,Y,Z,v,'--k','LineWidth', 0.8);
step = 1/norm;
v=-1*step:step:1*step;
v=[ 0 0];
[Co,h] = contour(X,Y,Z,v,'k','LineWidth', 0.8);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep'))
colormap cool


plot(x_(C1), y_(C1) , '.', 'Markersize', 10, 'Color', [1 0 0]);
plot(x_(C2), y_(C2) , '.', 'Markersize', 10, 'Color', 0.5*[1 1 1] );
%[sx,I] = sort(x_);
%plot(sx',fx(I)');
box    = PP*[min(x_) max(x_) min(y_) max(y_)];
axis(box);
set(gca,'FontSize',12);

%plot Support Vectors
plot(x_(SV1), y_(SV1) , 'o', 'Markersize', 6, 'Color', [1 0 0]);
plot(x_(SV2), y_(SV2) , 'o', 'Markersize', 6, 'Color', 0.5*[1 1 1]);



title(['C= ', num2str(C,4), ', ', 'Number of SVs: ', num2str(length(SV),3)]);

hold off;
drawnow;
toc;
