function [grid, P] = tauchen(N, m, rho, mu, sigma)
    %% grid
    thetaN = m * sigma / sqrt(1-rho^2);
    grid = linspace(-thetaN, thetaN, N);

    %% matriz de transicao
    delta = (max(grid) - min(grid)) / (N-1);
    [I,J] = meshgrid(grid, grid);
    P = normcdf((I + delta/2 - (1-rho)*mu - rho*J)/sigma) - ...
        normcdf((I - delta/2 - (1-rho)*mu - rho*J)/sigma);
    P(:,1) = normcdf((grid(1) - (1-rho)*mu - rho*grid + delta/2)/sigma);
    P(:,N) = 1 - normcdf((grid(N) - (1-rho)*mu - rho*grid - delta/2)/sigma);
    
end