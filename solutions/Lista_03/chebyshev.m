function [T,raizes] = chebyshev(K,j)
    T = cos(j*acos(K));
    i = (1:j)';
    raizes = -cos((2*i-1)/(2*j)*pi);
end
