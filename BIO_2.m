clear; clc; close all; 

%% punto 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretizzazione problema di Poisson con differenza finite al variare
% dell'ampiezza della mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 2.^(2:10); 
tspan = [0, 1];
alpha = 1; beta = -1; % dati al bordo

% soluzione esatta se f = 1

% mi aspetto che l'errore sia nulla dal momento che le differenze finite
% centrate per la derivata seconda sono esatte per polinomi di grado 3.

% f = @(x) 1;
% u = @(x) -0.5*x.^2 - 1.5*x + 1;

% soluzione esatta se f(x) = sin(x)
u = @(x) sin(x) - (2+sin(1))*x +1;
f = @(x) sin(x);

err = zeros(length(N), 1);
H = zeros(length(N), 1);
for j = 1:length(N)
    n = N(j);
    h = 1/(n+1);
    H(j) = h; 
    x = linspace(tspan(1), tspan(2), n+2).';
    
    A = -(1/h^2) * spdiags(repmat([1, -2, 1], n, 1), -1:1, n, n);
    
    F = f(x(2:end-1));
    F(1) = F(1) + alpha/h^2;
    F(n) = F(n) + beta/h^2;
    
    U = [alpha; A\F; beta];
    
    err(j) = max(abs(u(x)-U));
end

figure; 
loglog(H, err, LineWidth=2, DisplayName='$\|u-u_h\|_\infty$');
hold on
grid on
loglog(H, 1e-2*H.^2, '--r', DisplayName='$O(h^2)$')
legend(Location="best", Interpreter="latex")
title('Problema di Poisson con differenze finite', Interpreter='latex')


%% punto 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equazione del calore con differenze finite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

u_0t = 0; u_1t = 0;
u_x0 = @(x) 100*sin(pi*x);

u = @(t, x) 100 * exp(-(pi^2)*t) .* sin(pi*x); % soluzione esatta

T = 0.2;
m = 2^10; % tempo
n = 2^10; % spazio

U = zeros(n+2, m+1);

t = linspace(0, T, m+1).';
x = linspace(0, 1, n+2).';

U(:, 1) = u_x0(x);
U(1, :) = u_0t * ones(1, m+1);
U(n+2, :) = u_1t * ones(1, m+1);

Dt = T/m;
h = 1/(n+1);

A = -(1/h^2) * spdiags(repmat([1, -2, 1], n, 1), -1:1, n, n);
A = eye(n) + Dt*A;
dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

err = zeros(1, m);
for j=2:m+1
    
    F = U(2:n+1, j-1); % in questo problema f = 0
    
    U(2:end-1, j) = dA\F;
    err(j-1) = max(abs( u((j-1)*Dt, x) - U(:, j) ));

end

figure;
imagesc(t, x, U); colorbar eastoutside
title('Equazione del calore', Interpreter='latex')
xlabel('Intervallo temporale')
ylabel('intervallo spaziale')

