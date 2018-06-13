function [xn,a,b] = normalize(x)
    a = 2/range(x);
    b = - (max(x)+min(x)) / (max(x)-min(x));
    xn = a*x + b;
end