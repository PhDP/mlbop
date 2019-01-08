


function subgrad_coef = compute_subgrad_coef(x,d,e,n,i,epsilon,loss_type,loss_params)

switch loss_type
    case 'l2'
            if e(i)^2 <= epsilon 
                subgrad_coef = 0;
            else
                subgrad_coef = -2*e(i);
            end;
    case 'l1'
        if abs(e(i)) <= epsilon 
                subgrad_coef = 0;
            else
                subgrad_coef = -sign(e(i));
            end;
    case 'huber'
        sigma = loss_params(1);
        if abs(e(i)) <= epsilon 
                subgrad_coef = 0;
        else
            if abs(e(i)) <= sigma
                subgrad_coef = -e(i)/sigma;
            else
                subgrad_coef = -sign(e(i));
            end;
        end;
end;