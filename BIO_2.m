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
loglog(H, err, '-o', LineWidth=2, DisplayName='$\|u-u_h\|_\infty$');
hold on
grid on
loglog(H, err(length(N))*(H./H(length(N))).^2, '--r', DisplayName='$O(h^2)$')
legend(Location="best", Interpreter="latex")
title('Problema di Poisson con differenze finite', Interpreter='latex')


%% punto 2 (a)
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
imagesc(t, x, U); colorbar eastoutside; axis xy
title('Equazione del calore', Interpreter='latex')
xlabel('Intervallo temporale')
ylabel('intervallo spaziale')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Studio come l'errore dipende dal passo della discretizzazione Dt e h.
% In generale ho:
% errore del primo ordine in tempo; 
% errore del secondo ordine in spazio.
% Quindi si ha che errore = O(Dt) + O(h^2).
% Per studiare variazione errore in funzione di h impongo Dt = h^2
% Per studiare variazione errore in funzione di Dt impongo h = sqrt(Dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Studio errore in funzione di h (errore in spazio) %%%%%%%%%%%%%

N = 2.^(4:10);

H = zeros(length(N), 1);
err_h = zeros(length(N), 1);

for j = 1:length(N)
    n = N(j);

    h = 1/(n+1);
    Dt = h^2; % per studiare variazione errore in funzione di h
    m = ceil(T/Dt);

    H(j) = h;


    U = zeros(n+2, m+1);

    t = linspace(0, T, m+1).';
    x = linspace(0, 1, n+2).';

    U(:, 1) = u_x0(x);
    U(1, :) = u_0t * ones(1, m+1);
    U(n+2, :) = u_1t * ones(1, m+1);

    A = -(1/h^2) * spdiags(repmat([1, -2, 1], n, 1), -1:1, n, n);
    A = eye(n) + Dt*A;
    dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

    err = zeros(1, m);
    for jj=2:m+1

        F = U(2:n+1, jj-1); % in questo problema f = 0

        U(2:end-1, jj) = dA\F;
        err(jj-1) = max(abs( u((jj-1)*Dt, x) - U(:, jj) ));

    end

    err_h(j) = max(err);
    
end

figure; 
loglog(H, err_h, '-o', LineWidth=2, DisplayName='errore')
hold on
grid on
loglog(H, err_h(length(N))*(H./H(length(N))).^2, '--r', DisplayName='$O(h^2)$')
legend(Location="best", Interpreter='latex')
xlabel('h')
ylabel('errore')
title("Errore in funzione dell'ampiezza h", Interpreter='latex')


%%%%%%%%%%% Studio errore in funzione di Dt (errore in tempo) %%%%%%%%%%%%%

M = 2.^(4:10);
DT = zeros(length(M), 1);
err_t = zeros(length(M), 1);

for j = 1:length(M)
    m = M(j);
    Dt = T/m;
    DT(j) = Dt;
    h = sqrt(Dt); % per studiare variazione errore in funzione di Dt

    n = ceil(1/h)-1; 
    h = 1/(n+1);

    U = zeros(n+2, m+1);

    t = linspace(0, T, m+1).';
    x = linspace(0, 1, n+2).';

    U(:, 1) = u_x0(x);
    U(1, :) = u_0t * ones(1, m+1);
    U(n+2, :) = u_1t * ones(1, m+1);

    A = -(1/h^2) * spdiags(repmat([1, -2, 1], n, 1), -1:1, n, n);
    A = eye(n) + Dt*A;
    dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

    err = zeros(1, m);
    for jj=2:m+1

        F = U(2:n+1, jj-1); % in questo problema f = 0

        U(2:end-1, jj) = dA\F;
        err(jj-1) = max(abs( u((jj-1)*Dt, x) - U(:, jj) ));

    end

    err_t(j) = max(err);

end

figure; 
loglog(DT, err_t, '-o', LineWidth=2, DisplayName='errore')
hold on
grid on
loglog(DT, err_t(length(M))*DT./DT(length(M)), '--r', DisplayName='$O(\Delta t)$')
legend(Location="best", Interpreter="latex")
xlabel('$\Delta t$', Interpreter='latex')
ylabel('errore')

title("Errore in funzione dell'ampiezza $\Delta t$", Interpreter='latex')







