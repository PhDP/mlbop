

function value = kappa(x,y, kernel_type, kernel_params)

if (strcmp(kernel_type,'gaus')==1)
    sigma = kernel_params(1);
    N=length(x);
    norm = sum( (x-y).^2 );
    value = exp( -norm/sigma^2 );
elseif (strcmp(kernel_type,'gaus_c')==1)
    sigma = kernel_params(1);
    N=length(x);
    exponent = sum( (x-conj(y)).^2 );
    %value = 2*(exp( -real(exponent)/sigma^2 ));
    value = 2*real(exp(-exponent/sigma^2 ));
elseif (strcmp(kernel_type,'linear')==1)
    value = 0;
    N=length(x);
    for i=1:N
        value = value + x(i)*conj(y(i));
    end;
elseif (strcmp(kernel_type,'poly')==1)
    d = kernel_params(1);
    value = (1 + x*transpose(y))^d;
    %value = ( (1 + x*transpose(y))/( sqrt(real(x*transpose(x)) * real(y*transpose(y)) ) ) )^d;
elseif (strcmp(kernel_type,'poly_c')==1)
    d = kernel_params(1);
    value = 2*real( (1 + x*transpose(y))^d );
    %value = 2*real( ( (1 + x*transpose(y))/( sqrt(real(x*transpose(x)) * real(y*transpose(y)) ) ) )^d );
end;
    
