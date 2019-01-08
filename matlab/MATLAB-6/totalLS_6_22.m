%==========================================================
% Exercise 6.22
% This function generates the model and the data and 
% runs the Total Least Squares 
% and the Least Squares algorithms.
%==========================================================


N=150; %Number of observations
l=90;  %problem size
sv=0.01; % Variance of the additive noise
sv2=0.2;
Iter=100;
 
sx=size(sv);
t1=zeros(l,1);
t2=zeros(l,1);
t3=zeros(l,1);
t4=zeros(l,1);
theta=randn(l,1);
%%%%%%%%%%%%%%%%INPUT NOISE FREE SCENARIO
for i=1:Iter
A=randn(N,l);

noise=sqrt(sv)*randn(N,1);
y=A*theta;

y1=y+noise;
A1=A;
that = inv(A1'*A1)*A1'*y1;
t1=t1+that;
 
%%%%%%%%%%%%ADD SOME NOISE ON THE INPUT%%%%%%%%%%
E=sv2*randn(N,l);
A2=A+E;

that2 = inv(A2'*A2)*A2'*y1;%Compute the classical Least Squares
 
t2=t2+that2;

%%%%%%%%TOTAL LEAST SQUARES%%%%%%%%%%
Ahatext=[A2 y1];
[U S V]=svd(Ahatext);%Compute the SVD of the Matrix
sl1= min(S(find(S~=0)));%Find the smallest eigenvalue not equal to zero.
that3 = inv(A2'*A2-sl1^2*eye(l))*A2'*y1;
t3=t3+that3;

Ahatext=[A2 y1];
that4 = tls(A2,y1,eps);
t4=t4+that4;

end
norm(t1./Iter-theta)
norm(t2./Iter-theta)
norm(t3./Iter-theta)
norm(t4./Iter-theta)