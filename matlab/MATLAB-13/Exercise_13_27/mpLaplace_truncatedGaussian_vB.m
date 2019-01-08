function [w_vb, bt]= mpLaplace_truncatedGaussian_vB(Phi, y, Phi_gram, MaxIter, const)
% Function 
% Variational Bayes for linear regression based on a multiparameter 
% truncated Laplace prior
% 
% Input
% Phi:      Endmember matrix
% y:        Pixel's spectral measurements
% Phi_gram: Gramian matrix of Phi, Phi.' * Phi 
% MaxIter:  Maximum Iterations ~ 1e3
% const:    small constant ~ 1e-8
% 
% Output
% w_vb:     Estimated abundance vector
% bt:       Estimated noise precision
% 
% Copyright K. Themelis, A. Rontogiannis, K. Koutroumbas

[M,N] = size(Phi);

% indexes
indx = zeros(N-1,N); for n = 1 : N, tmp = 1:N; tmp(n) = []; indx(:,n) = tmp; end

alph = ones(N,1); % initialization
bn = ones(N,1);
mu = zeros(N,1);
mutr = zeros(N,1);
sgmitr = ones(N,1);

t = 0;
while 1
    t = t + 1;
    mutr_old = mutr;

    bt = (2e-6 + N + M ) / (2e-6 + sum(alph .* (mutr.^2 + sgmitr)) + norm(y-Phi*mutr) + sgmitr.' * diag(Phi_gram) );

    % untruncated variance
    sgmi =  (alph + diag(Phi_gram) ).^-1 / bt;
    sgmisd = sqrt(sgmi);
    sgmitr = sgmi;
    
    for n = 1 : N
        % untruncated mean
        mu(n) =  Phi(:,n).' * (y - Phi(:,indx(:,n)) * mutr(indx(:,n))) / (alph(n) + Phi_gram(n,n));
        
        % truncated mean
        tmp = 1 - .5 * erfc(mu(n) /sqrt(2) /sgmisd(n) );
        if tmp > const
            mutr(n) = mu(n) + sgmisd(n) / sqrt(2 * pi) * exp(-.5 * mu(n)^2 / sgmi(n) ) / tmp; % truncated mean
            tmp1 = exp( - mu(n)^2 / sgmi(n) / 2)  / sqrt(2*pi)  / tmp ;
            sgmitr(n) = sgmi(n) * (1 - mu(n) / sgmisd(n) * tmp1 - tmp1^2 ); % truncated variance
        else
            if mu(n) < 0
                mu(n) = 1e-6;
            end
        end
    end
        
    alph = sqrt( bn ./ (bt * (mutr.^2 + sgmitr)));
    
    bn = (2e-6 + 2) ./ (2e-6 + 1./alph + 1./bn );
    
    if ( norm(mutr - mutr_old) < const || t > MaxIter), w_vb = mu; break; end
end