%*******************************************************************
%SMO algorithm for regression.
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
%kernel_type:    The type of the kernel, usually this ['gaus', 'poly']
%kernel_params   The parameters for the kernel (i.e. sigma, e.t.c.)
%-----------------------------------------------------------------
%
%-----------------------------------------------------------------
%output variables
%
% a  :          a vector of size m. These are the support vectors
% 
% b  :          a real number. This is the offset of the support vector
%                   expansion.



function [a1, a2, b, KKT, rmse] = SMO_regression2(x, y, C, epsilon, kernel_type, kernel_params)

global  a1_g a2_g b_g u_g KKT_g NB_g;
tol=0.001;

[M,N]=size(x);


if strcmp(kernel_type, 'gaus')
    par = kernel_params;
    norms=zeros(M,M);
    for i=1:M
        T = bsxfun(@minus,x,x(i,:));
        norms(i,:) = sum(T.^2,2);
    end;
    Kernel_matrix(:,:) = exp(-norms./(par^2));
elseif strcmp(kernel_type, 'gaus_c')
    par = kernel_params;
    norms=zeros(M,M);
    for i=1:M
        T = bsxfun(@minus,x,conj(x(i,:)));
        norms(i,:) = sum(T.^2,2);
    end;
    Kernel_matrix(:,:) = 2*real(exp(-norms./(par^2)));
else
    for i=1:M
        for j=1:M
            Kernel_matrix(i,j) = kappa(x(i,:), x(j,:), kernel_type, kernel_params);
        end;
    end;
end;

    
    
%initialize the support vectors
a1_g = zeros(M,1);
a2_g = zeros(M,1);


%initialize the threshold
b_g=0;

%u contains the values of the sv expansion for each input
u_g = b_g*ones(M,1);

%KKT(i) is 1 if the i-th sv (a(i)) satisfies the KKT conditions
%KKT(i) is 0 otherwise
KKT_g=zeros(M,1);
%Update KKT conditions
update_KKT(x,y,C,epsilon);


%NB(i) is 1 if the sv a(i) is non bound, i.e. 0<a(i)<C
%NB(i) is 0 otherwise.
NB_g=zeros(M,1);

numChanged = 0;
examineAll = 1;
LoopCounter = 0;
while ( numChanged > 0) || (examineAll==1)
    LoopCounter = LoopCounter+1;
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
    if (mod(LoopCounter,2)==0)
        MinimumNumChanged = max(1, 0.1*M);
    else
        MinimumNumChanged=1;
    end;
    if (examineAll == 1)
        examineAll = 0;
    else if (numChanged < MinimumNumChanged)
            examineAll = 1;
        end;
    end;
    
end;
a1=a1_g;
a2=a2_g;
b=b_g;
KKT=KKT_g;

rmse = sqrt( mean( (u_g - y).^2 ) );








function ret=examineExample(j, x, y, C, epsilon, Kernel_matrix)

global  KKT_g NB_g;
M=length(y);


if ( KKT_g(j) == 0 )
    if ( sum(NB_g)>1 )
        %find the i - second choice heurestic
        i = second_choice_heurestic(j,x,y,epsilon);
        if ( takeStep(i, j,  x, y, C, epsilon, 1, Kernel_matrix) )
            ret = 1;
            return;
        end;
    end;
    
    %loop over all non-bound a's starting at a random point
    i0 = randi(M);
    for i = i0 : i0+M-1
        if NB_g(mod(i-1,M)+1)==1
            if (takeStep(mod(i-1,M)+1, j,  x, y, C, epsilon, 0, Kernel_matrix))
                ret=1;
                return;
            end;
        end;
    end;
    
    %loop over all possible a's starting at a random point
    i0 = randi(M);
    for i = i0 : i0+M-1
        if (takeStep(mod(i-1,M)+1, j,  x, y, C, epsilon, 0, Kernel_matrix)==1)
            ret=1;
            return;
        end;
    end;
end;
ret = 0;
    
    



function ret = takeStep(i,j, x, y, C, epsilon, chosen, Kernel_matrix)
global  a1_g a2_g b_g u_g NB_g KKT_g;
global  a1i_new a1j_new a2i_new a2j_new;

tol = 0.001;

if (i == j)
    ret=0;
    return;
end;
E_i =  y(i) - u_g(i);
E_j =  y(j) - u_g(j);
gamma = a1_g(i) + a1_g(j) - a2_g(i) - a2_g(j);


kii = Kernel_matrix(i,i);
%kii = kappa(x(i,:), x(i,:), kernel_type, kernel_params);
kij = Kernel_matrix(i,j);
%kij = kappa(x(i,:), x(j,:), kernel_type, kernel_params);
kjj = Kernel_matrix(j,j);
%kjj = kappa(x(j,:), x(j,:), kernel_type, kernel_params);
eta = -2*kij + kii + kjj;

case1=0;
case2=0;
case3=0;
case4=0;
finished=0;
a1i_old=a1_g(i);
a2i_old=a2_g(i);
a1j_old=a1_g(j);
a2j_old=a2_g(j);

while (finished==0)
    if (case1==0) && ( (a1_g(i)>0) || ((a2_g(i)==0)&&(E_i-E_j>0)) ) && ( (a1_g(j)>0)||( (a2_g(j)==0)&&(E_i-E_j)<0) )
        L = max(0, gamma-C);
        H = min(gamma,C);
        if (L<H)
            a2 = a1_g(j) - (E_i-E_j)/eta;
            a2 = min(a2,H);
            a2=max(L,a2);
            a1 = a1_g(i) - (a2 - a1_g(j));
            if abs(a1_g(i) - a1) + abs(a1_g(j) - a2) > tol
                a1_g(i) = a1;
                a1_g(j) = a2;
            end;
        else
            finished = 1;
        end;
        case1 = 1;
    elseif (case2==0) && ( (a1_g(i)>0) || ((a2_g(i)==0) && (E_i-E_j>2*epsilon)) ) && ( (a2_g(j)>0)||((a1_g(j)==0)&&(E_i-E_j>2*epsilon)) )
        L=max(0,-gamma);
        H=min(C,C-gamma);
        if (L<H)
            a2 = a2_g(j) + (E_i-E_j-2*epsilon)/eta;
            a2 = min(a2,H);
            a2 = max(L,a2);
            a1 = a1_g(i) + (a2-a2_g(j));
            if abs(a1_g(i) - a1) + abs(a2_g(j) - a2) > tol
                a1_g(i) = a1;
                a2_g(j) = a2;
            end;
        else
            finished = 1;
        end;
        case2=1;
    elseif (case3 == 0) && ((a2_g(i)>0)||((a1_g(i)==0)&&(E_i-E_j<-2*epsilon))) && ((a1_g(j)>0)||((a2_g(j)==0)&&(E_i-E_j<-2*epsilon)))
        L=max(0,gamma);
        H=min(gamma+C,C);
        if (L<H)
            a2 = a1_g(j) - (E_i-E_j+2*epsilon)/eta;
            a2=min(a2,H);
            a2=max(L,a2);
            a1 = a2_g(i) + (a2 - a1_g(j));
            if abs(a2_g(i) - a1) + abs(a1_g(j) - a2) > tol
                a2_g(i) = a1;
                a1_g(j) = a2;
            end;
        else
            finished = 1;
        end;
        case3=1;
    elseif (case4==0) && ( (a2_g(i)>0)||((a1_g(i)==0)&&(E_i-E_j)<0) ) && ( (a2_g(j)>0)||((a1_g(j)==0)&&(E_i-E_j>0)) )
        L=max(0,-gamma-C);
        H=min(-gamma,C);
        if (L<H)
            a2 = a2_g(j) + (E_i-E_j)/eta;
            a2 = min(a2,H);
            a2 = max(L,a2);
            a1 = a2_g(i) - (a2 - a2_g(j));
            if abs(a2_g(i) - a1) + abs(a2_g(j) - a2) > tol
                a2_g(i) = a1;
                a2_g(j) = a2;
            end;
        else
            finished = 1;
        end;
        case4=1;
    else
        finished = 1;
    end;

    %update errors
    E_i = E_i - ((a1_g(i) - a2_g(i)) - (a1i_old - a2i_old))*Kernel_matrix(i,i) - ((a1_g(j) - a2_g(j)) - (a1j_old - a2j_old))*Kernel_matrix(i,j);
    E_j = E_j - ((a1_g(i) - a2_g(i)) - (a1i_old - a2i_old))*Kernel_matrix(j,i) - ((a1_g(j) - a2_g(j)) - (a1j_old - a2j_old))*Kernel_matrix(j,j); 
end;


if abs(a1_g(i)-a1i_old) + abs(a1_g(j)-a1j_old) + abs(a2_g(i)-a2i_old) + abs(a2_g(j)-a2j_old) < tol
    ret = 0;
    return;
end;
ret = 1;

%update b
M=length(y);
b_old=b_g;
b1=0;
b2=0;
N1=0;
N2=0;
for k=1:M
    if (a1_g(k)>0) && (a1_g(k)<C)
        %b = y_i - <w,x_i> - e
        in_prod = 0;
        for l=1:M
            in_prod = in_prod + (a1_g(l) - a2_g(l))*Kernel_matrix(l,k);
        end;
        b1 = b1  +  y(k) - in_prod - epsilon;
        N1=N1+1;
        %break;
    end;
    
    if (a2_g(k)>0) && (a2_g(k)<C)
        %b = y_i - <w,x_i> + e
        in_prod = 0;
        for l=1:M
            in_prod = in_prod + (a1_g(l) - a2_g(l))*Kernel_matrix(l,k);
        end;
        b2 = b2 + y(k) - in_prod + epsilon;
        N2=N2+1;
        %break;
    end;
end;

if (N1>0)&&(N2>0)
    b_g = (b1/N1 + b2/N2)/2;
elseif N1>0
    b_g = b1/N1;
elseif N2>0
    b_g = b2/N2;
else
    b1=0;
    in_prod = 0;
    for l=1:M
        in_prod = in_prod + (a1_g(l) - a2_g(l))*Kernel_matrix(l,k);
    end;
    b1 = b1  +  y(k) - in_prod;
    N1=N1+1;
    b_g = b1/N1;
end;
        

%Update u vector
M=length(y);
for k=1:M
    %u_g(k) = 0;
    u_g(k) = u_g(k) - b_old + b_g + (a1_g(i) - a2_g(i) - a1i_old + a2i_old)*Kernel_matrix(k,i) + (a1_g(j) - a2_g(j) - a1j_old + a2j_old)*Kernel_matrix(k,j);
    %for l=1:M
    %    u_g(k) = u_g(k) + (a1_g(l) - a2_g(l))*kappa(x(k,:), x(l,:));
    %end;
    %u_g(k) = u_g(k) + b_g;
end;


%Update KKT conditions
update_KKT(x, y, C, epsilon)

%Update NB vector
for k=1:M
    if ( (a1_g(k)>0) && (a1_g(k)<C) ) || ( (a2_g(k)>0) && (a2_g(k)<C) )
       NB_g(k)=1; %Non Bound example
    else
        NB_g(k)=0; %Bound example
    end;
end;














function  i = second_choice_heurestic(j,x,y,epsilon)
global   u_g;
M=length(y);

E_1 = y(1) - u_g(1);
E_j = y(j) - u_g(j);
max1 = abs(  E_1 - E_j  );
i1=1;
for k=2:M
    E_k = y(k) - u_g(k);
    if abs( E_k - E_j ) > max1
        max1 = abs( E_k - E_j );
        i1=k;
    end;
end;

max2 = abs( E_j - E_1 + 2*epsilon );
i2=1;
for k=2:M
    E_k = y(k) - u_g(k);
    if abs( E_j - E_k + 2*epsilon ) > max2
        max2 = abs( E_j - E_k + 2*epsilon );
        i2=k;
    end;
end;

max3 = abs( E_1 - E_j + 2*epsilon );
i3=1;
for k=2:M
    E_k = y(k) - u_g(k);
    if abs( E_k - E_j + 2*epsilon ) > max3
        max3 = abs( E_k - E_j + 2*epsilon );
        i3=k;
    end;
end;

max4 = abs( E_j - E_1 );
i4=1;
for k=2:M
    E_k = y(k) - u_g(k);
    if abs( E_j - E_k ) > max4
        max4 = abs( E_j - E_k );
        i4=k;
    end;
end;

i=i1;
maximum=max1;
if max2>maximum
    maximum = max2;
    i=i2;
end;
if max3>maximum
    maximum=max3;
    i=i3;
end;
if max4>maximum
    maximum=max4;
    i=i4;
end;


        


function update_KKT(x, y, C, epsilon)
global  a1_g a2_g u_g KKT_g;

M=length(y);
tol = 0.01;

for i=1:M
    E_i = u_g(i) - y(i);
    if ( (E_i>=epsilon) && (a2_g(i)<C) ) || ( (E_i<epsilon)&&(a2_g(i)>0) ) || ( (-E_i>=epsilon)&&(a1_g(i)<C) ) || ( (-E_i<epsilon)&&(a1_g(i)>0) )
        KKT_g(i)=0;
    else
        KKT_g(i)=1;
    end;

 end;






