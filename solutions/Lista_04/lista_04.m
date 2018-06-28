% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 1.0001;
rho = 0.9;
sigma = 0.01;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);

%% a) 
%  Use the Tauchen method to discretize the stochastic process in a
%  Markov chain with 9 states. (Use 3 standard deviations for each side)

[zgrid, P] = tauchen(9,3,rho,0,sigma);
zgrid;
P

%% b)
% Discretize the asset space using a grid and solve the individual goat 
% farmer problem for each state variable.

egrid = exp(zgrid) % endowment grid

r = 1/q - 1;
min_a = 0;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000)'; % asset grid

[E, A, AP] = meshgrid(egrid, agrid, agrid); % maximizing over 3rd dim

% Consumo e utilidade
C = E + A - q*AP;
C(C<0) = 0;
U = u(C);

% Força bruta
V0 = repmat(agrid, 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [3,2,1]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 3);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
subplot(4,2,1)
plot(agrid, V)
title('Função valor')
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,2)
plot(agrid, agrid(idx))
hold on
plot([0 60], [0 60], '--k')
hold off
title('Função política')
xlabel('a')
ylabel('G(a,z)') 

%% c)
% Find the stationary distribution Pi(z, a) and use it to compute the 
% aggregate savings in the economy.

disp('c) -------------------------------------------------------------');

Pi0 = repmat(1/9000, 1000, 9);

I = zeros(1000, 9, 1000);
for i = 1:1000
    I(:,:, i) = (agrid(idx) == agrid(i));
end

Pi = zeros(1000,9);
err = 1;
while err > 10^-5
    for i = 1:1000
        Pi(i,:) = sum((Pi0 .* I(:,:,i))*P);
    end
    err = norm(Pi-Pi0);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));

fprintf('Aggregate savings: %2.2f\n', S);


%% d)
% Suppose rho = .97. Redo the analysis. How does the savings rate 
% compare now? Explain.
disp('d) -------------------------------------------------------------');
% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 1.0001;
rho = 0.97;
sigma = 0.01;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);

[zgrid, P] = tauchen(9,3,rho,0,sigma);
P

egrid = exp(zgrid) % endowment grid

r = 1/q - 1;
min_a = 0;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000)'; % asset grid

[E, A, AP] = meshgrid(egrid, agrid, agrid); % maximizing over 3rd dim

% Consumo e utilidade
C = E + A - q*AP;
C(C<0) = 0;
U = u(C);

% Força bruta
V0 = repmat(agrid, 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [3,2,1]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 3);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

Pi0 = repmat(1/9000, 1000, 9);

I = zeros(1000, 9, 1000);
for i = 1:1000
    I(:,:, i) = (agrid(idx) == agrid(i));
end

Pi = zeros(1000,9);
err = 1;
while err > 10^-5
    for i = 1:1000
        Pi(i,:) = sum((Pi0 .* I(:,:,i))*P);
    end
    err = norm(Pi-Pi0);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));

fprintf('Aggregate savings: %2.2f\n', S);

% Gráficos
subplot(4,2,3)
plot(agrid, V)
title('Função valor')
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,4)
plot(agrid, agrid(idx))
hold on
plot([0 60], [0 60], '--k')
hold off
title('Função política')
xlabel('a')
ylabel('G(a,z)') 


%% e)
% Suppose gamma = 5. Redo the analysis. How does the savings rate compare now? 
% Explain.

disp('e) -------------------------------------------------------------');
% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 5;
rho = 0.9;
sigma = 0.01;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);

[zgrid, P] = tauchen(9,3,rho,0,sigma);
P

egrid = exp(zgrid) % endowment grid

r = 1/q - 1;
min_a = 0;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000)'; % asset grid

[E, A, AP] = meshgrid(egrid, agrid, agrid); % maximizing over 3rd dim

% Consumo e utilidade
C = E + A - q*AP;
C(C<0) = 0;
U = u(C);

% Força bruta
V0 = repmat(agrid, 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [3,2,1]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 3);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

Pi0 = repmat(1/9000, 1000, 9);

I = zeros(1000, 9, 1000);
for i = 1:1000
    I(:,:, i) = (agrid(idx) == agrid(i));
end

Pi = zeros(1000,9);
err = 1;
while err > 10^-5
    for i = 1:1000
        Pi(i,:) = sum((Pi0 .* I(:,:,i))*P);
    end
    err = norm(Pi-Pi0);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));

fprintf('Aggregate savings: %2.2f\n', S);

% Gráficos
subplot(4,2,5)
plot(agrid, V)
title('Função valor')
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,6)
plot(agrid, agrid(idx))
hold on
plot([0 60], [0 60], '--k')
hold off
title('Função política')
xlabel('a')
ylabel('G(a,z)') 



%% f)
% Suppose sigma = .05. Redo the analysis. How does the savings rate compare 
% now? Explain.
disp('f) -------------------------------------------------------------');
% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 1.0001;
rho = 0.9;
sigma = 0.05;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);

[zgrid, P] = tauchen(9,3,rho,0,sigma);
P

egrid = exp(zgrid) % endowment grid

r = 1/q - 1;
min_a = 0;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000)'; % asset grid

[E, A, AP] = meshgrid(egrid, agrid, agrid); % maximizing over 3rd dim

% Consumo e utilidade
C = E + A - q*AP;
C(C<0) = 0;
U = u(C);

% Força bruta
V0 = repmat(agrid, 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [3,2,1]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 3);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

Pi0 = repmat(1/9000, 1000, 9);

I = zeros(1000, 9, 1000);
for i = 1:1000
    I(:,:, i) = (agrid(idx) == agrid(i));
end

Pi = zeros(1000,9);
err = 1;
while err > 10^-5
    for i = 1:1000
        Pi(i,:) = sum((Pi0 .* I(:,:,i))*P);
    end
    err = norm(Pi-Pi0);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));

fprintf('Aggregate savings: %2.2f\n', S);

% Gráficos
subplot(4,2,7)
plot(agrid, V)
title('Função valor')
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,8)
plot(agrid, agrid(idx))
hold on
plot([0 60], [0 60], '--k')
hold off
title('Função política')
xlabel('a')
ylabel('G(a,z)') 
