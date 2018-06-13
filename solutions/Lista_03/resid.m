function [R] = resid(gamma, K0, alpha, beta, delta, mu)
    C0 = C_hat(gamma, K0);
    K1 = K0.^alpha + (1-delta)*K0 - C0;
    C1 = C_hat(gamma, K1);
    R = beta * (C1./C0).^(-mu) .* (1-delta+alpha*K1.^(alpha-1));    
end


function [C] = C_hat(gamma, K)
    % gamma vetor coluna (d,1)
    % K vetor coluna (n,1)
    j = size(gamma,1);
    T = zeros(size(K,1), j);
    for i = 1:j
        T(:,i) = chebyshev(K,i);
    end
    C = T * gamma;
end