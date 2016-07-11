%*******************************************************************
%SMO algorithm for classification.
%started: 22/10/2010
%ended: 22/10/2010
%-------------------------------------------------------------------
%Description
%
%u_i = Ö(x_i) = Ó y_n a_n K(x_i,x_n) + b
%
% The Algorithm computes the parameters a_k, of the expansion of the
% solution w = Sum a_n*K(. , x_n) and the parameter b.
%-------------------------------------------------------------------
%input variables
%
%x:              a M x N matrix. It contains the input vectors x_i (each one in
%                   every row.
%y:              a vector of size M. The class of each input vector 
%C:              the trade-off parameter of the SVM model
%-----------------------------------------------------------------
%
%-----------------------------------------------------------------
%output variables
%
% a  :          a vector of size m. These are the support vectors
% 
% b  :          a real number. This is the offset of the support vector
%                   expansion.



function [a, b] = SMO_classification(x, y, C, epsilon, kernel_type, kernel_params)

global  a_g b_g u_g KKT_g NB_g;
tol=0.001;

[M,N]=size(x);


if kernel_type == 'gaus'
    par = kernel_params;
    norms=zeros(M,M);
    for i=1:M
        T = bsxfun(@minus,x,x(i,:));
        norms(i,:) = sum(T.^2,2);
    end;
    Kernel_matrix(:,:) = exp(-norms./(par^2));
else
    for i=1:M
        for j=1:M
            Kernel_matrix(i,j) = kappa(x(i,:), x(j,:), kernel_type, kernel_params);
        end;
    end;
end;



%initialize the support vectors
a_g = zeros(M,1);

%initialize the threshold
b_g=0;

%u contains the values of the sv expansion for each input
u_g = zeros(M,1);

%KKT(i) is 1 if the i-th sv (a(i)) satisfies the KKT conditions
%KKT(i) is 0 otherwise
KKT_g=zeros(M,1);
%Update KKT conditions
for k=1:M
    r2 = (u_g(k)-y(k))*y(k);
    if ( (r2<-tol)&&(a_g(k)<C) ) || ( (r2>tol) && (a_g(k)>0) )
        KKT_g(k)=0; %KKT condition not satisfied
    else
        KKT_g(k)=1;
    end;
end;

%NB(i) is 1 if the sv a(i) is non bound, i.e. 0<a(i)<C
%NB(i) is 0 otherwise.
NB_g=zeros(M,1);

numChanged=0;
examineAll=1;
while ( numChanged > 0) || (examineAll==1)
    numChanged=0;
    if (examineAll==1)
        %loop over all training examples
        for i=1:M
            numChanged = numChanged + examineExample(i, x, y, C, epsilon, Kernel_matrix);
        end;
    else
        %loop over all training examples, where a(i) is non-bound
        for i=1:M
            if (NB_g(i) == 1)
                numChanged = numChanged + examineExample(i, x, y, C, epsilon, Kernel_matrix);
            end;
        end;
    end;
    if (examineAll == 1)
        examineAll = 0;
    else if (numChanged == 0)
            examineAll = 1;
        end;
    end;
    
end;
a=a_g;
b=b_g;
disp('SMO Finished');









function ret=examineExample(j, x, y, C, epsilon, Kernel_matrix)

global  u_g KKT_g NB_g;
M=length(y);

E_j = u_g(j) - y(j);
if ( KKT_g(j) == 0 )
    if ( sum(NB_g)>1 )
        %find the i - second choice heurestic
        i = second_choice_heurestic(j,x,y);
        if ( takeStep(i, j,  x, y, C, epsilon, Kernel_matrix) )
            ret = 1;
            return;
        end;
    end;
    
    %loop over all non-bound a's starting at a random point
    i0 = randi(M);
    for i = i0 : i0+M-1
        if NB_g(mod(i-1,M)+1)==1
            if (takeStep(mod(i-1,M)+1, j,  x, y, C, epsilon, Kernel_matrix))
                ret=1;
                return;
            end;
        end;
    end;
    
    %loop over all possible a's starting at a random point
    i0 = randi(M);
    for i = i0 : i0+M-1
        if (takeStep(mod(i-1,M)+1, j,  x, y, C, epsilon, Kernel_matrix)==1)
            ret=1;
            return;
        end;
    end;
end;
ret = 0;
    
    







function ret = takeStep(i,j, x, y, C, epsilon, Kernel_matrix)
global  a_g b_g u_g KKT_g NB_g;

tol = 0.001;

if (i == j)
    ret=0;
    return;
end;
E_i = u_g(i) - y(i);
E_j = u_g(j) - y(j);
s = y(i)*y(j);
%compute L, H
if ( y(i)~=y(j) )
    L = max(0, a_g(j)-a_g(i) );
    H = min(C, C+a_g(j)-a_g(i));
else
    L = max(0, a_g(j)+a_g(i)-C );
    H = min(C, a_g(j)+a_g(i));
end;
if (L==H)
    ret = 0;
    return;
end;

kii = Kernel_matrix(i,i);
kij = Kernel_matrix(i,j);
kjj = Kernel_matrix(j,j);
eta = kii + kjj - 2*kij;

if (eta > 0)
    aj_new = a_g(j) + y(j)*(E_i - E_j)/eta;
    if (aj_new < L)
        aj_new = L;
    else if (aj_new>H)
            aj_new = H;
        end;
    end;
else
    %when numerical errors are involved
    kij = Kernel_matrix(i,j);
    f1 = y(i)*(E_i+b_g) - a_g(i)*Kernel_matrix(i,i) - s*a_g(j)*Kernel_matrix(i,j);
    f2 = y(j)*(E_j+b_g) - s*a_g(i)*Kernel_matrix(i,j) - a_g(j)*Kernel_matrix(j,j);
    L1 = a_g(i) + s*(a_g(j) - L);
    H1 = a_g(i) + s*(a_g(j) - H);
    Psi_L = L1*f1 + L*f2 + 0.5*L1^2*Kernel_matrix(i,i) + 0.5*L^2*Kernel_matrix(j,j) + s*L*L1*Kernel_matrix(i,j);
    Psi_H = H1*f1 + H*f2 + 0.5*H1^2*Kernel_matrix(i,i) + 0.5*H^2*Kernel_matrix(j,j) + s*H*H1*Kernel_matrix(i,j);
    L_obj = Psi_L;
    H_obj = Psi_H;
    if (L_obj < H_obj - epsilon)
        aj_new = L;
    else
        if (L_obj > H_obj + epsilon)
            aj_new = H;
        else
            aj_new = a_g(j);
        end;
    end;
end;

if ( abs(aj_new - a_g(j)) < epsilon*(aj_new + a_g(j) + epsilon) )
    ret = 0;
    return;
end;

ai_new = a_g(i) + s*( a_g(j) - aj_new);



%Update threshold b
b1 = E_i + y(i)*(ai_new - a_g(i))*Kernel_matrix(i,i) + y(j)*(aj_new - a_g(j))*Kernel_matrix(i,j) + b_g;
b2 = E_j + y(i)*(ai_new - a_g(i))*Kernel_matrix(i,j) + y(j)*(aj_new - a_g(j))*Kernel_matrix(j,j) + b_g;
if (ai_new > 0)&&(ai_new < C)
    b_g = b1;
elseif (aj_new > 0)&&(aj_new < C)
    b_g = b2;
else
    b_g = (b1+b2)/2;
end;

ret=1;
%Update a vector
a_g(i) = ai_new;
a_g(j) = aj_new;


%Update u vector
M=length(y);
for k=1:M
    u_g(k) = 0;
    for l=1:M
        u_g(k) = u_g(k) + y(l)*a_g(l)*Kernel_matrix(k,l);
    end;
    u_g(k) = u_g(k) - b_g;
end;


%Update KKT conditions
for k=1:M
    r2 = (u_g(k)-y(k))*y(k);
    if ( (r2<-tol)&&(a_g(k)<C) ) || ( (r2>tol) && (a_g(k)>0) )
        KKT_g(k)=0; %KKT condition not satisfied
    else
        KKT_g(k)=1;
    end;
end;

%Update NB vector
for k=1:M
    if (a_g(k)>0) && (a_g(k)<C)
       NB_g(k)=0; %Bound example
    else
        NB_g(k)=1; %non Bound example
    end;
end;




function  i = second_choice_heurestic(j,x,y)
global  a_g b_g u_g KKT_g NB_g;

maximum = abs( u_g(1) - y(1) - ( u_g(j) - y(j) ));
i=1;
M=length(y);
for k=1:M
    if abs( u_g(k) - y(k) - ( u_g(j) - y(j) ) ) > maximum
        maximum = abs( u_g(k) - y(k) - ( u_g(j) - y(j) ) );
        i=k;
    end;
end;
        

    

