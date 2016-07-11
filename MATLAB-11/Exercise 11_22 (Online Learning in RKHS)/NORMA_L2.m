%*******************************************************************
%NORMA  algorithm with L2 norm.
%started: 21/07/2014
%ended: 21/07/2014
%updated: 21/07/2014
%-------------------------------------------------------------------
%Description
%
%y_n = <Phi(u_n), w>
%minimize   |d_n - y_n|^2 + ë||f||^2
%
% The Algorithm computes the parameters a_k, of the expansion of the
% solution w = Sum a_k*K(z_n, .)
% Employs Quantization
%-------------------------------------------------------------------
%input variables
%
%z:              a N x L matrix. z_k = (x_{k-L+1}, x_{k-L+2}, . . ., x_{k})
%d:              a N x 1 matrix. 
%N:              an integer number.
%mu:             a real number (the step size).
%print_flag:     1 if we want to print the MSE 
%                0 if we don't want to print the MSE
%kernel_type:    'gaus', for the gaussian kernel
%                'poly', for the polynomial kernel (d=2)
%                'vovk_inf'
%                'vovk_poly'
%                'diric', for the dirichlet kernel
%                'period'
%                 SEE ALSO kappa.m

%kernel_params:     a vector of real numbers. It contains all the required
%                   parameters for the kernel evaluation. For example the sigma for the
%                   gaussian kernel
%-----------------------------------------------------------------
%
%-----------------------------------------------------------------
%output variables
%
% a                 :   a matrix containing the real expansion coefficients
% centers           :   a N0 x L matrix, containing
%                       the z_n that contribute to the expansion
% e                 :   a vector containing the error at each instance
% expansion_size    :   Nx1 vector containing the number of the centers at
%                       each time instant. It can be used to visualize  
%                       how the size of the expansion evolves.



function [a, centers, e, expansion_size] = QNORMA_L2(z,d,N,mu,lambda,print_flag,kernel_type,kernel_params,sparse_params)

epsilon=0.000001;

%----------------------------------------------------------------
%quantization size
q_size = sparse_params(1);  %typically 1
%--------------------------------------------------------------


weight = 1;
n0 = sparse_params(1);

N0=0;%the number of the centers that appear in the expansion
for n=1:N
    %compute y_n and e_n
    %---------------------------------------------------
    %for slow matalb computations
    if strcmp(kernel_type,'gaus_f')==0
        R_sum=0;
        for k=1:N0
            R_sum = R_sum + a(k)*Kappa(centers(k,:)',z(n,:)',kernel_type,kernel_params);
        end;
    else
    %------------------------------------------------
    %to speed up things
    %we use a C-based procedure to compute the filter output
    %this works only with gaussian kernel
        if (N0==0)
            R_sum=0;
        else
            T = bsxfun(@minus, centers, z(n,:));
            norms = sum(T.^2,2);
            Kernel_matrix = exp(-norms./(kernel_params(1)^2));
            R_sum = a*Kernel_matrix;
        end;
    end;
    
     %------------------------------------------------
    y(n) = R_sum ;
    e(n) = d(n) - y(n);
     
    %sparsification rules
    %sparsification rules
    
    
    add_this_center=1;
    
    %find the minimum distance between the new center z(n,:) and all
    %the previous ones.
    if n<=n0
        add_only = 1;
    else
        add_only = 0;
    end;
        
    
    if (add_only == 1)
        N0=N0+1;
        centers(N0,:) = z(n,:);
        a(N0) = +2*mu*e(n);
        a(1:N0-1) = (1-lambda*mu)*a(1:N0-1);
    else
        centers(1:n0-1,:) = centers(2:n0,:);
        centers(N0,:) = z(n,:);
        a(1:n0-1) = (1-lambda*mu)*a(2:n0);
        a(n0) = +2*mu*e(n);
    end;
        
    expansion_size(n) = N0;
    
    mu = mu*weight;
    
end;

plot_v = e;
if print_flag==1
    plot(10*log10(sum_mean_val(abs(plot_v(length(z(1,:))+1:N)).^2)),'-r','LineWidth',3);
end;
    
    
        
    

