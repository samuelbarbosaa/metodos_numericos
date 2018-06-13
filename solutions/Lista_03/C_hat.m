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
        