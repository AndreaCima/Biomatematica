clear; clc; close all; 

Cm = 1; 
GNa = 120;
GK = 36;
GL = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;
sigma = 1e-3; 

alpha_m = @(V) 0.1*(25-V).*(exp((25-V)/10)-1).^(-1);
alpha_h = @(V) 0.07*exp(-V/20);
alpha_n = @(V) 0.01*(10-V) .* (exp( (10-V)/10 ) -1 ).^(-1);

beta_m = @(V) 4*exp(-V/18);
beta_h = @(V) (exp( (30-V)/10 ) + 1).^(-1);
beta_n = @(V) 0.125*exp(-V/80);

V0 = @(x) 6.6*zeros(length(x), 1);
m0 = @(x) 5.2934e-2*zeros(length(x), 1);
h0 = @(x) 5.9611e-1*zeros(length(x), 1);
n0 = @(x) 3.1768e-1*zeros(length(x), 1);
Iapp =@(t, x) 20 * (0 <= t & t <= 1 & 0 <= x & x <= 0.05);

T = 20;

n_t = 1000; % tempo
n_x = 1000; % spazio

Dt = T/n_t;
step_h = 1/(n_x+1);

x = linspace(0, 1, n_x+2)';
t = linspace(0, T, n_t+1)';

V = zeros(n_x+2, n_t+1); V(:, 1) = V0(x);
m = zeros(n_x+2, n_t+1); m(:, 1) = m0(x);
h = zeros(n_x+2, n_t+1); h(:, 1) = h0(x);
n = zeros(n_x+2, n_t+1); n(:, 1) = n0(x);


A_V = -(1/step_h^2) * spdiags(repmat([1, -2, 1], n_x+2, 1), -1:1, n_x+2, n_x+2);
A_V(1, 2) = -2/(step_h^2); A_V(n_x+2, n_x+1) = -2/(step_h^2); % valori di correzione per matrice A
A_V = (Cm/Dt).*speye(n_x+2) + sigma*A_V; 
dA_V = decomposition(A_V); % per risolvere il sistema lineare in modo più efficiente


for j = 2:n_t+1
    F_V = Cm/Dt * V(:, j-1) - ( GNa * m(:, j-1).^3 .* h(:, j-1) .* (V(:, j-1)-VNa) + ...
        GK * n(:, j-1).^4 .* (V(:, j-1) - VK) + GL*(V(:, j-1)-VL)) + Iapp(t(j), x);
    V(:, j) = dA_V \ F_V;

    m(:, j) = (m(:, j-1) / Dt + alpha_m(V(:, j-1))) ./ (1/Dt + alpha_m(V(:, j-1)) + beta_m(V(:, j-1)));
    h(:, j) = (h(:, j-1) / Dt + alpha_h(V(:, j-1))) ./ (1/Dt + alpha_h(V(:, j-1)) + beta_h(V(:, j-1)));
    n(:, j) = (n(:, j-1) / Dt + alpha_n(V(:, j-1))) ./ (1/Dt + alpha_n(V(:, j-1)) + beta_n(V(:, j-1)));
end

figure; 
mesh(t, x, V)

