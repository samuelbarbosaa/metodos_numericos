clear;

%% Ambiente
beta   = 0.987;
mu     = 2;
alpha  = 1/3;
delta  = 0.012;
rho    = 0.95;
sigma  = 0.007;

%% Estado estacionario
kss = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1))

%% Grids
[zgrid, P] = tauchen(7, 3, rho, 0, sigma);
zgrid = exp(zgrid);
P

kgrid = linspace(0.75*kss, 1.25*kss, 500)';
[kgridn, a, b] = normalize(kgrid);
plot(kgridn)

%% Computing C1
d = 4;
gamma0 = ones(d+1,1);
[~, K0] = chebyshev(0,d+1);
R = @(gamma) resid(gamma, K0, kgrid, alpha, beta, delta, mu, a, b);
gamma_otimo = fsolve(R, gamma0);
C = C_hat(gamma_otimo, kgridn);

K1 = kgrid.^alpha + (1-delta)*kgrid - C;

plot(kgrid,K1)