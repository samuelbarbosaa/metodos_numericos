function [R] = resid(gamma, K0, Z, P, alpha, beta, delta, mu, a, b)
    C0 = C_hat(gamma, K0);
    K = (K0-b)/a;
    K1 = Z.*K.^alpha + (1-delta)*K - C0;
    K1N = a*K1 + b;
    C1 = C_hat(gamma, K1N);
    R = beta * (C1./C0).^(-mu) .* (1-delta+alpha*Z.*K1.^(alpha-1)) * P' - 1;    
end


