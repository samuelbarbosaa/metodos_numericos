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
[z_grid, P] = tauchen(7, 3, rho, 0, sigma);
z_grid = exp(z_grid)
P

k_grid = linspace(0.75*kss, 1.25*kss, 500);
kp_grid = k_grid;
[K, KP, Z] = meshgrid(k_grid,kp_grid,z_grid);

%% Consumo e utilidade
Y = Z .* K .^ alpha;
S = KP - (1-delta).*K;
C = Y - S;
C(C<0) = 0.001;
U = (C.^(1-mu)-1)/(1-mu);

% Chute inicial para a função valor em cada um dos 7 estados
V0 = repmat(kp_grid', 1, 7);

%% Força bruta
it=1;
err=1;
tol=10^-5;
itmax=1000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,500), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = max(max(abs(V-V0)));
    if mod(it,20)==0 disp(err); end
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
plot(k_grid, V)

figure()
plot(kp_grid, kp_grid(idx))


%% Monotonicidade
V0 = repmat(kp_grid', 1, 7);

% Primeira iteração
EV = V0*P';
EV = permute(repmat(EV,1,1,500), [1,3,2]);
H = U + beta * EV;
[V, idx0] = max(H, [], 1);
V=squeeze(V); idx0=squeeze(idx0);
V0=V;

it=1;
err=1;
tol=10^-5;
itmax=1000;

% Iterações monótonas
tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,500), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H(idx0:end,:,:), [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = max(max(abs(V-V0)));
    if mod(it,20)==0 disp(err); end
    it=it+1;
    V0=V;
    if idx0 <= 500 idx0=idx+1; end
end
toc()

% Gráficos
figure()
plot(k_grid, V)

figure()
plot(kp_grid, kp_grid(idx))


%% Concavidade
V0 = repmat(kp_grid', 1, 7);

it=1;
err=1;
tol=10^-5;
itmax=1000;
check=1;

tic()
while check && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,500), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    check = any(any(V>=V0));
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
plot(k_grid, V)

figure()
plot(kp_grid, kp_grid(idx))

%% Acelerador



%% Multigrid



%% Grid endógeno
