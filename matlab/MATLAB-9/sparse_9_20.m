clc; clear all

A=[0.5 2 1.5;1.5 2.3 3.5];
x=[2.5 0 0]';

y = A*x;

% L2 norm minimized solution
theta2=A'*inv(A*A')*y;
error_L2 = norm(y-A*theta2); %Check that theta2 is a solution
disp('L2 norm minimization solution:')
disp(theta2) % Note that thata2 is not one-sparse

disp(['Error achieved with L2 norm minimization:',num2str(error_L2)])
disp('------------------------')
% Exhaustive search of all combinations
disp('Start checking for potential 1-sparse solutions')
disp('Check solution [x,0,0]')
subA = A(:,1);
xx1=zeros(3,1);
xx1(1)=inv(subA'*subA)*subA'*y;
disp('Result:')
disp(xx1)
error1 = norm(y-A*xx1); %Check that xx1 is a solution
disp(['Achieved error:', num2str(error1)])

disp('------------------------')
disp('Check solution [0,x,0]')
subA = A(:,2);
xx2=zeros(3,1);
xx2(2)=inv(subA'*subA)*subA'*y;
disp('Result:')
disp(xx2)
error2 = norm(y-A*xx2); %Check that xx2 is a solution
disp(['Achieved error:', num2str(error2)])

disp('------------------------')
disp('Check solution [0,0,x]')
subA = A(:,3);
xx3=zeros(3,1);
xx3(3)=inv(subA'*subA)*subA'*y;
disp('Result:')
disp(xx3)
error3 = norm(y-A*xx3); %Check that theta2 is a solution
disp(['Achieved error:', num2str(error3)])

disp('------------------------')
disp('So, the ector ')
disp(xx1)
disp('is a solution so no reason to search for 2-sparse solutions')
theta0 = xx1;

disp('------------------------')
disp('Check the L2 norms of the L2 and L0 minimization solutions')
disp(['L2 minimization: ',num2str(norm(theta2))])
disp(['L0 minimization: ',num2str(norm(theta0))])

% the theta2 has smaller L2 norm than theta0 as expected.



