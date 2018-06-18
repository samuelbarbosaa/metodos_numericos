function [C] = C_hat(gamma, K)
    % gamma matriz (d,z)
    % K matriz (n,z)
    j = size(gamma,1);
    T = zeros(size(K,1), j, size(K,2));
    for i = 1:j
        T(:,i,:) = chebyshev(K,i);
    end
    % T(k,j,z)
    
    for z = 1:size(K,2)
        C(:,z) = T(:,:,z) * gamma;
    end
end
        