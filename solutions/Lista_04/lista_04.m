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
zgrid
P

%% b)
% Discretize the asset space using a grid and solve the individual goat 
% farmer problem for each state variable.

egrid = exp(zgrid); % endowment grid

r = 1/q - 1;
min_a = - egrid(1) / r;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000); % asset grid

[AP, A, E] = meshgrid(agrid, agrid, egrid);

% Consumo e utilidade
C = E + A - q*AP;
C(C<=0) = 0.001;
U = u(C);

% Força bruta
V0 = repmat(agrid', 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

% Gráficos
figure()
subplot(4,2,1)
plot(agrid, V(1000:-1:1,:))
title('Função valor')
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,2)
plot(agrid, agrid(idx))
hold on
plot([-24,54], [-24,54], 'k--')
hold off
title('Função política')
xlabel('a')
ylabel('G(a,z)') 

%% c)
% Find the stationary distribution Pi(z, a) and use it to compute the 
% aggregate savings in the economy.

id0 = repmat(1/9, 1, 9);
err = 1;
while err > 10^-5
    id = id0 * P;
    err = norm(id-id0);
    id0 = id;
end


Pi0 = repmat(1/(1000*9), 1000, 9);

err = 1;
while err > 10^-5
    Pi  = Pi0 * P;
    err = norm(Pi0 - Pi);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));
fprintf('c) Aggregate savings: %2.2f\n', S);


%% d)
% Suppose rho = .97. Redo the analysis. How does the savings rate 
% compare now? Explain.

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

egrid = exp(zgrid); % endowment grid

r = 1/q - 1;
min_a = - egrid(1) / r;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000); % asset grid

[AP, A, E] = meshgrid(agrid, agrid, egrid);

C = E + A - q*AP;
C(C<=0) = 0.001;
U = u(C);

V0 = repmat(agrid', 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

subplot(4,2,3)
plot(agrid, V(1000:-1:1,:))
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,4)
plot(agrid, agrid(idx))
hold on
plot([-22, 55], [-22,55], 'k--')
hold off
xlabel('a')
ylabel('G(a,z)') 

Pi0 = repmat(1/(1000*9), 1000, 9);

err = 1;
while err > 10^-5
    Pi  = Pi0 * P;
    err = norm(Pi0 - Pi);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));
fprintf('d) Aggregate savings: %2.2f\n', S);



%% e)
% Suppose gamma = 5. Redo the analysis. How does the savings rate compare now? 
% Explain.


% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 5;
rho = 0.97;
sigma = 0.01;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);
[zgrid, P] = tauchen(9,3,rho,0,sigma);

egrid = exp(zgrid); % endowment grid

r = 1/q - 1;
min_a = - egrid(1) / r;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000); % asset grid

[AP, A, E] = meshgrid(agrid, agrid, egrid);

C = E + A - q*AP;
C(C<=0) = 0.001;
U = u(C);

V0 = repmat(agrid', 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

subplot(4,2,5)
plot(agrid, V(1000:-1:1,:))
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,6)
plot(agrid, agrid(idx))
hold on
plot([-22,55], [-22,55], 'k--')
hold off
xlabel('a')
ylabel('G(a,z)') 

Pi0 = repmat(1/(1000*9), 1000, 9);

err = 1;
while err > 10^-5
    Pi  = Pi0 * P;
    err = norm(Pi0 - Pi);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));
fprintf('e) Aggregate savings: %2.2f\n', S);



%% f)
% Suppose sigma = .05. Redo the analysis. How does the savings rate compare 
% now? Explain.

% Ambiente
clear

beta = 0.96;
q = 0.96;
gamma = 5;
rho = 0.97;
sigma = 0.05;

u = @(c) (c.^(1-gamma)-1)/(1-gamma);
up = @(c) c.^(-gamma);
[zgrid, P] = tauchen(9,3,rho,0,sigma);

egrid = exp(zgrid); % endowment grid

r = 1/q - 1;
min_a = - egrid(1) / r;
max_a = 2 * egrid(9) / r;

agrid = linspace(min_a, max_a, 1000); % asset grid

[AP, A, E] = meshgrid(agrid, agrid, egrid);

C = E + A - q*AP;
C(C<=0) = 0.001;
U = u(C);

V0 = repmat(agrid', 1, 9);

it=1;
err=1;
tol=10^-5;
itmax=3000;

tic()
while err>tol && it<itmax
    EV = V0*P';
    EV = permute(repmat(EV,1,1,1000), [1,3,2]);
    H = U + beta * EV;
    [V, idx] = max(H, [], 1);
    V=squeeze(V); idx=squeeze(idx);
    err = norm(V-V0);
    it=it+1;
    V0=V;
end
toc()

subplot(4,2,7)
plot(agrid, V(1000:-1:1,:))
xlabel('a')
ylabel('V(a,z)') 

subplot(4,2,8)
plot(agrid, agrid(idx))
hold on
plot([-13, 89], [-13,89], 'k--')
hold off
xlabel('a')
ylabel('G(a,z)') 

Pi0 = repmat(1/(1000*9), 1000, 9);

err = 1;
while err > 10^-5
    Pi  = Pi0 * P;
    err = norm(Pi0 - Pi);
    Pi0 = Pi;
end

S = sum(sum(Pi .* agrid(idx)));
fprintf('f) Aggregate savings: %2.2f\n', S);
