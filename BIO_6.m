% Modello FHN in 1D: propagazione del potenziale d'azione lungo una fibra

clear; clc; close all; 

sigma = 1e-3;
b = 5;
c = 1;
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

v0 = @(x) zeros(length(x), 1);
w0 = @(x) zeros(length(x), 1);

m = 2^10;
n = 2^10;

T = 50; 

Dt = T/m;
h = 1/(n+1);

t = linspace(0, T, m+1)';
x = linspace(0, 1, n+2)';

V = zeros(n+2, m+1);
V(:, 1) = v0(x);

W = zeros(n+2, m+1);
W(:, 1) = w0(x);

A = -(1/h^2) * spdiags(repmat([1, -2, 1], n+2, 1), -1:1, n+2, n+2);
A(1, 2) = -2/(h^2); A(n+2, n+1) = -2/(h^2); % valori di correzione per matrice A
A = (1/Dt).*speye(n+2) + sigma*A; 
dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

I_app = @(t, x) 2*(0 <= x & x <= 0.04 & 0 <= t & t <= 1);


for j=2:m+1
    F = V(:, j-1)/Dt + b * V(:, j-1) .* ( V(:, j-1)-beta ) .* ( delta - V(:, j-1) ) - c*W(:, j-1) + I_app(t(j), x);

    V(:, j) = dA\F;
    W(:, j) = (Dt / (1+Dt*e*gamma)) * (W(:, j-1)./Dt + e*V(:, j-1));
end

figure;
imagesc(t, x, V); axis xy
xlabel('tempo')
ylabel('spazio')
title('Potenziale', Interpreter='latex')

figure;
mesh(t, x, V);
xlabel('tempo')
ylabel('spazio')
title('Potenziale', Interpreter='latex')


figure;
imagesc(t, x, W); axis xy
xlabel('tempo')
ylabel('spazio')
title('Recovery', Interpreter='latex')


figure;
mesh(t, x, W); axis xy
xlabel('tempo')
ylabel('spazio')
title('Recovery', Interpreter='latex')










