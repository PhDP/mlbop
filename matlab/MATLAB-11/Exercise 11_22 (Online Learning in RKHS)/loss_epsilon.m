


function val = loss_epsilon(ksi, epsilon, loss_type, loss_params)

switch loss_type
    case 'l2'
        val = max(0, ksi^2 - epsilon);
    case 'l1'
        val = max(0, abs(ksi) - epsilon);
    case 'huber'
        sigma = loss_params(1);
        if abs(ksi) <= sigma
            val = max(0, ksi^2/(2*sigma) - epsilon);
        else
            val = max(0, abs(ksi) - sigma/2 - epsilon);
        end;
end;