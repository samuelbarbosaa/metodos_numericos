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
    err = norm(V-V0);
    if mod(it,20)==0 disp(err); end
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
subplot(1,2,1)
plot(k_grid, V)
title('Função valor')
xlabel('k')
ylabel('V(k,z)') 

subplot(1,2,2)
plot(kp_grid, kp_grid(idx))
title('Função política')
xlabel('k')
ylabel('G(k,z)') 

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
    err = norm(V-V0);
    if mod(it,20)==0 disp(err); end
    it=it+1;
    V0=V;
    if idx0 <= 500 idx0=idx; end
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
    check = any(any(V<V0));
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
it=1;
err=1;
tol=10^-5;
itmax=1000;

V0 = repmat(kp_grid', 1, 7);

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,500), [1,3,2]);
    H = U + beta * EV;
    if it < 20 || mod(it, 50)==0  
        [V, idx] = max(H, [], 1);
        V=squeeze(V);idx=squeeze(idx);
        err = max(max(abs(V-V0)));
    else
        for j=1:7
            for i=1:500
                V(i,j)=H(idx(i,j),i,j);
            end
        end
        err = max(max(abs(V-V0)));
    end
    if mod(it,50)==0 disp(err); end
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
plot(k_grid, V)

figure()
plot(kp_grid, kp_grid(squeeze(idx)))



%% Multigrid

%% 100 pontos
k100 = linspace(0.75*kss, 1.25*kss, 100);
kp100 = k100;
[K, KP, Z] = meshgrid(k100,kp100,z_grid);

% Consumo e utilidade
Y = Z .* K .^ alpha;
S = KP - (1-delta).*K;
C = Y - S;
C(C<0) = 0.001;

u = @(c,mu) (c.^(1-mu)-1)/(1-mu);
U = u(C,mu);

it=1;
err=1;
tol=10^-5;
itmax=1000;

V0 = repmat(kp100', 1, 7);

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,100), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = max(max(abs(V-V0)));
    if mod(it,20)==0 disp(err); end
    it=it+1;
    V0=V;
end
toc()

%% 500 pontos

k500 = linspace(0.75*kss, 1.25*kss, 500);
kp500 = k500;
[K, KP, Z] = meshgrid(k500,kp500,z_grid);

% Consumo e utilidade
Y = Z .* K .^ alpha;
S = KP - (1-delta).*K;
C = Y - S;
C(C<0) = 0.001;

u = @(c,mu) (c.^(1-mu)-1)/(1-mu);
U = u(C,mu);

it=1;
err=1;
tol=10^-5;
itmax=1000;

V0 = interp1(k100, V, k500, 'linear');

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

%% 5000 pontos

k5000 = linspace(0.75*kss, 1.25*kss, 5000);
kp5000 = k5000;
[K, KP, Z] = meshgrid(k5000,kp5000,z_grid);

% Consumo e utilidade
Y = Z .* K .^ alpha;
S = KP - (1-delta).*K;
C = Y - S;
C(C<0) = 0.001;

u = @(c,mu) (c.^(1-mu)-1)/(1-mu);
U = u(C,mu);

it=1;
err=1;
tol=10^-5;
itmax=1000;

V0 = interp1(k500, V, k5000, 'linear');

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,5000), [1,3,2]);
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
plot(k5000, V)

figure()
plot(kp5000, kp5000(squeeze(idx)))


%% Grid endógeno

clear;

% Ambiente
beta   = 0.987;
mu     = 2;
alpha  = 1/3;
delta  = 0.012;
rho    = 0.95;
sigma  = 0.007;

% Estado estacionario
kss = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1))

% Grids
[z_grid, P] = tauchen(7, 3, rho, 0, sigma);
z_grid = exp(z_grid)
P

k_grid = linspace(0.75*kss, 1.25*kss, 500);
kp_grid = k_grid;
[K, KP, Z] = meshgrid(k_grid,kp_grid,z_grid);

% Consumo e utilidade
Y = Z .* K .^ alpha;
S = KP - (1-delta).*K;
C = Y - S;
C(C<0) = 0.001;

u = @(c,mu) (c.^(1-mu)-1)/(1-mu);
U = u(C,mu);

% Finding optimal k' from k
G0 = squeeze(repmat(1:500, 1, 1, 7));

for j=1:7
    for i=1:500
        KPP=H(idx(i,j),i,j);
    end
end


C.^(-mu) .* alpha .* Z .* KP.^(alpha-1) + 1 - delta
