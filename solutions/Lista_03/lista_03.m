clear;

%% Ambiente
beta   = 0.987;
mu     = 2;
alpha  = 1/3;
delta  = 0.012;
rho    = 0.95;
sigma  = 0.007;

%% Estado estacionario
kss = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));

%% Grids
[zgrid, P] = tauchen(7, 3, rho, 0, sigma);
zgrid = exp(zgrid);
kgrid = linspace(0.75*kss, 1.25*kss, 500)';
[kgridn, a, b] = normalize(kgrid);
[Z,K] = meshgrid(zgrid, kgrid);
[~,KN] = meshgrid(zgrid, kgridn);

%% Computing G(k,z) and C(k,z);

% d = 4
% gamma0 = [3 1 0 0 0]';

d = 6;

gamma0 = rand(d+1,1);

[~, cg] = chebyshev(0,d+1);
[Z0, K0] = meshgrid(zgrid, cg);
R = @(gamma) resid(gamma, K0, Z0, P, alpha, beta, delta, mu, a, b);

opts = optimoptions('fsolve', 'MaxFunctionEvaluations', 5000);
gamma_otimo = fsolve(R, gamma0, opts);

C = C_hat(gamma_otimo, KN); 
K1 = Z.*K.^alpha + (1-delta)*K - C; % Politica capital

CP = C;
for i=1:7
    CP(:,i) = interp1(K(:,i), C(:,i), K1(:,i)); % Politica consumo
end

subplot(1,3,1)
plot(kgrid,CP)
title 'Função política do consumo - C(k,z)'

subplot(1,3,2)
plot(kgrid,K1)
title 'Função política de capital - G(k,z)'


%% Euler errors
up = @(c) c.^(-mu);
up_inv = @(c) c.^(-1/mu);
pmgK = @(K, Z) alpha*Z.*K.^(alpha-1) + 1 - delta;

E = up(CP).* pmgK(K1, Z)*P';
EE = log10(abs(1 - up_inv(beta*E)./C));

subplot(1,3,3)
plot(kgrid,EE)
title 'Euler errors'
