



function value=Kappa(x,y,type,param)

switch type
    case 'gaus'
        sigma=param(1);
        N=length(x);
        norm=0;
        for i=1:N
            norm = norm + (x(i)-y(i))^2;
        end;
        value = exp( -norm/sigma^2 );
        
     case 'gaus_f'
        sigma=param(1);
        N=length(x);
        norm=0;
        for i=1:N
            norm = norm + (x(i)-y(i))^2;
        end;
        value = exp( -norm/sigma^2 );
    
    case 'poly'
        d=param(1);
        value=(transpose(x)*y + 1)^d;
        
    case 'vovk_inf'
        p=param(1);
        value = (1 - transpose(x)*y)^p / (1 - transpose(x)*y);
        
    case 'vovk_poly'
        value = 1 / (1 - transpose(x)*y);
        
    case 'diric'
        n=param(1);
        N=length(x);
        norm=0;
        for i=1:N
            norm = norm + (x(i)-y(i))^2;
        end;
        value = sin((2*n+1)*norm/2)/(sin(norm/2)+0.0001);
        
    case 'period'
        n=param(1);
        sigma=param(2);
        N=length(x);
        norm=0;
        for i=1:N
            norm = norm + (x(i)-y(i))^2;
        end;
        sum=0;
        for k=0:n
            sum = sum+exp(-k^2*sigma^2/2)*cos(k*norm);
        end;
        value=sum;
end;