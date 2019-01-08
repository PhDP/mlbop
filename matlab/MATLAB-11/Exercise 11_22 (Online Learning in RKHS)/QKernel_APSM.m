%*******************************************************************
%Kernel APSM algorithm with quantization sparsification method.
%started: 05/02/2011
%ended: 11/02/2011
%Created by Pantelis Bouboulis
%-------------------------------------------------------------------
%Description
%
%Phi(u_n) ------> filter ----> dhat_n = < Phi(x_n), w> 
%
% The Algorithm computes the parameters a_k, of the expansion of the
% solution w = Sum a_k*K(z_n, .) 
%-------------------------------------------------------------------
%input variables
%
%x:         a N x L matrix. x_k = (x_{k-L+1}, x_{k-L+2}, . . ., x_{k})
%d:         a N x 1 matrix. 
%N:         an integer number (the number of input samples).
%epsilon:   the value for the epsilon-defined cost.
%Q:         the number of samples concurrently processed.

%loss_type:      'l2' for the l2 cost function
%                'l1' for the l1 cost function
%                'huber' for the huber loss function
%loss_params:    the parameters for the loss function, if needed
%                for example the sigma for the huber loss.
%kernel_type:    'gaus', for the gaussian kernel
%                'gaus_f' for the gaussian kernel with fast C coding 
%                'poly', for the polynomial kernel
%                'vovk_inf'
%                'vovk_poly'
%                'diric', for the dirichlet kernel
%                'period'
%                 SEE ALSO Kappa.m
%kernel_params:      a vector of real numbers. It contains all the required
%                   parameters for the kernel evaluation. For example the sigma for the
%                   gaussian kernel


%sparse_params:  a vector containing the parameters for the sparsification
%                 [Delta q_size]

%print_flag:     1 if we want to print the MSE 
%                0 if we don't want to print the MSE


%-----------------------------------------------------------------
%
%-----------------------------------------------------------------
%output variables
%
% a         :   a vector containing the  expansion coefficients
% e         :   a vector containing the error at each instance
%centers    :   a vector containing the expansion centers
%               for sparsification with the l2 ball (or no sparsification)
%               this is identical with  x



function [a,centers,e,expansion_size] = QKernel_APSM(x,d,N,epsilon,Q,loss_type,loss_params,kernel_type,kernel_params,sparse_params,print_flag)

a=[];
e=zeros(N,1);
expansion_size = zeros(N,1);
N0=0;
centers=[];

d_hat  = zeros(Q,1);
b = zeros(Q,1);
subgrad = zeros(Q,1);
norm_subgrad = zeros(Q,1);
omega = zeros(Q,1);
e_n = zeros(Q,1);
L_n = zeros(Q,1);
norm_f=0;

mu_n=0.5;
weight=1;
minimum_mu = 0.001;


Delta = sparse_params(1);
q_size = sparse_params(2);




pos = zeros(Q,1);
for n=1:N
    
    
    %find the minimum distance between the new center x(n,:) and all
    %the previous ones.
    minimum=1000000000;
    thesis=0;
    for k=1:N0
        dist = sum(abs(x(n,:)-centers(k,:)).^2);
        if (dist < minimum)
            minimum=dist;
            thesis = k;
        end;
    end;
    if (minimum < q_size) %do not add this center
        add_this_center=0;
    else
        add_this_center=1;
    end;
    
    
    for i=Q:-1:2
        pos(i) = pos(i-1);
    end;
    
    
    if (add_this_center==1)
        a(N0+1)=0;
        N0=N0+1;
        centers(N0,:) = x(n,:);
        pos(1) = N0;
    else
        pos(1) = thesis;
    end;
    expansion_size(n)=N0;
    
    
    
    q = min(N0,Q);
    for i=1:q
        %compute the (n-i+1)-th filter output
        %---------------------------------------------------
        %for slow matlab computations
        d_hat(i) = 0;
        if strcmp(kernel_type,'gaus_f')==0
            for j=1:n-1
                d_hat(i) = d_hat(i) + a(j)*Kappa(centers(j,:),x(n-i+1,:),kernel_type,kernel_params);
            end;
        %------------------------------------------------
        %to speed up things
        %we use a C-based procedure to compute the filter output
        %this works only with gaussian kernel
        else
            d_hat(i) = fast_real_output_kernel_computation(a, centers', x(n-i+1,:), N0, kernel_params(1));
        end;
        %------------------------------------------------
        
        e_n(i) = d(n-i+1) - d_hat(i);
        L_n(i) = loss_epsilon(e_n(i), epsilon, loss_type, loss_params);
        %compute the subgradient (coefficient) of the n-i+1 loss function
        subgrad(i) = compute_subgrad_coef(centers,d,e_n,n,i,epsilon,loss_type,loss_params);
        %compute the norm of each subgradient
        if subgrad(i)~=0
            %norm_subgrad(i) = subgrad(i)^2 * Kappa(x(n-i+1,:), x(n-i+1,:), kernel_type, kernel_params);
            norm_subgrad(i) = subgrad(i)^2 * fast_kernel(x(n-i+1,:), x(n-i+1,:), kernel_params(1));
        else
            norm_subgrad(i) = 0;
        end;   
    end;
    %store the error at instance n.
    e(n) = e_n(1);
    %how many subgradients are non zero?
    r = non_zero_cardinality(subgrad);
    %compute the omega coefficients
    if r>0
        for i=1:q
            omega(i)=1/r;
        end;
    else
        for i=1:q
            omega(i)=0;
        end;
    end;
    
    
    
    %update the current estimate
    %first coefficients
    %a(N0+1)=0;
    b = zeros(Q,1);
    for i=1:q
        if norm_subgrad(i) >= 0.000001
            b(q-i+1) = omega(i) *  L_n(i) * subgrad(i) / norm_subgrad(i);
            %a(N0+1-i+1) = a(N0+1-i+1) - mu_n * b(q-i+1);
        end;
    end;
    
    %extrapolation parameter
%     nominator = 0;
%      for i=1:q
%         if norm_subgrad(i) >= 0.000001
%             nominator = nominator + omega(i) *  L_n(i)*L_n(i) / norm_subgrad(i);
%         end;
%     end;
%     denominator=0;
%     for i=1:q
%         for j=1:q
%             denominator = denominator + b(q-i+1)*b(q-j+1)*Kappa(x(n-i+1,:), x(n-j+1,:), kernel_type, kernel_params);
%         end;
%     end;
%     if (denominator > 0.0001)
%         mu_n=nominator / denominator;
%     else
%         mu_n = 1;
%     end;
    
    
    mu_n = mu_n*weight;
    if mu_n<minimum_mu;
        mu_n=minimum_mu;
        weight=1;
    end;
    
    
    
    %update the current estimate
    for i=1:q
        if norm_subgrad(i) >= 0.000001
            a(pos(i)) = a(pos(i)) - mu_n * b(q-i+1);
        end;
    end;
    
   
   %L2-ball projection
    if mod(n,500)==0
        if strcmp(kernel_type,'gaus_f')==0
            norm_f = function_norm(a,centers,n,kernel_type,kernel_params);
        else
            norm_f = fast_compute_function_norm(a, centers', N0, kernel_params(1));
        end;

        if norm_f > Delta
            a = Delta/norm_f * a;
            norm_f=Delta;
        end;
    end;
   
        
        
end;

    
plot_v = e;
if print_flag==1
    plot(10*log10(sum_mean_val(abs(plot_v(length(x(1,:))+1:N)).^2)),'-g','LineWidth',3);
end;
