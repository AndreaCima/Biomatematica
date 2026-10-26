clear; clc; close all; 

%% Equazione di Kolmogorov - Fisher
% equazione di Kolmogorov-Fisher (diffusione-reazione) con condizioni di
% Neumann omogenee

%%%%%%%%%%%%%%%%%%%%%%%%%    PARAMETRI    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 1e-3;
b = 5;
alpha = 0.1;
 
I_app = 0;

T = 10;
m = 2^10; % tempo
n = 2^10; % spazio 

u_x0 = @(x) alpha;

t = linspace(0, T, m+1).';
x = linspace(0, 1, n+2).';

U = zeros(n+2, m+1);

U(:, 1) = u_x0(x); % condizione iniziale in spazio

% passi della discretizzazione 
Dt = T/m;
h = 1/(n+1);

A = -(1/h^2) * spdiags(repmat([1, -2, 1], n+2, 1), -1:1, n+2, n+2);

A(1, 2) = -2/(h^2); A(n+2, n+1) = -2/(h^2); % valori di correzione per matrice A

% matrice dei sistemi lineari che devo risolvere ad ogni istante temporale
A = (1/Dt).*speye(n+2) + sigma*A; 
dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

for j=2:m+1
    F = U(:, j-1)/Dt + b * U(:, j-1) .* ( 1 - U(:, j-1) ) + I_app;

    U(:, j) = dA\F;
end


figure; 

mesh(t, x, U);
xlabel('Intervallo temporale', Interpreter='latex')
ylabel('Intervallo spaziale', Interpreter='latex')

title('Kolmogorov - Fisher', Interpreter='latex')


% per alpha = 1 e per alpha = 0 (e I_app=0) ho soluzioni stazionaria dal
% momento che si annulla il termine di reazione e abbiamo visto che per la
% matrice di Neumann il vettore 1 è nel nucleo dal momento che la somma
% lungo le righe è nulla 


%%
clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% altra configurazione dell'equazione di Kolmogorov-Fisher, sempre con
% condizioni di Neumann omogenee.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_x0 = 0;
alpha = 10;
sigma = 1e-3;
b = 5;

I_app = @(t, x) alpha*(0 <= x & x <= 0.04 & 0 <= t & t <= 1);

T = 10;
m = 200; % tempo
n = 100; % spazio 


t = linspace(0, T, m+1).';
x = linspace(0, 1, n+2).';

U = zeros(n+2, m+1);

U(:, 1) = u_x0; % condizione iniziale in spazio

Dt = T/m;
h = 1/(n+1);

A = -(1/h^2) * spdiags(repmat([1, -2, 1], n+2, 1), -1:1, n+2, n+2);

A(1, 2) = -2/(h^2); A(n+2, n+1) = -2/(h^2); % valori di correzione per matrice A

% matrice dei sistemi lineari che devo risolvere ad ogni istante temporale
A = (1/Dt).*speye(n+2) + sigma*A; 
dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

for j=2:m+1
    F = U(:, j-1)/Dt + b * U(:, j-1) .* ( 1 - U(:, j-1) ) + I_app(t(j), x);

    U(:, j) = dA\F;
end


% animazione in tempo della soluzione sull'intervallo [0, 1]
figure; 

for j = 1:m+1
    plot(x, U(:, j))
    grid on
    ylim([0 2])
    title(sprintf('Tempo: t = %.2f s', t(j)), Interpreter="latex");
    xlabel('intervallo spaziale', Interpreter='latex')
    if j == 1
        pause(2) % la prima volta ci va un po' di tempo perchè apra il plot
    end
    pause(0.005)
end

make_video(U, [0 2])



