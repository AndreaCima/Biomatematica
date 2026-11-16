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

m = 800; % tempo
n = 1000; % spazio

T = 50; 

% passi temporali e spaziali
Dt = T/m;
h = 1/(n+1);

t = linspace(0, T, m+1)';
x = linspace(0, 1, n+2)';

V = zeros(n+2, m+1); % potenziale
V(:, 1) = v0(x);

W = zeros(n+2, m+1); % recovery
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
xlabel('tempo (ms)')
ylabel('spazio')
title('Potenziale', Interpreter='latex')

figure;
mesh(t, x, V);
xlabel('tempo (ms)')
ylabel('spazio')
title('Potenziale', Interpreter='latex')


figure;
imagesc(t, x, W); axis xy
xlabel('tempo (ms)')
ylabel('spazio')
title('Recovery', Interpreter='latex')


figure;
mesh(t, x, W); axis xy
xlabel('tempo (ms)')
ylabel('spazio')
title('Recovery', Interpreter='latex')
 

figure; % plot evoluzione temporale nel punto medio

y_min = min(min(V(:)), min(W(:))) - 0.05;
y_max = max(max(V(:)), max(W(:))) + 0.05;
idx_mid = ceil((n+2)/2); % indice punto medio

h1 = animatedline('LineWidth', 1.5, Color='blue', DisplayName='potenziale');
h2 = animatedline('LineWidth', 1.5, Color='red', DisplayName='recovery');

grid on
hold on
xlabel('tempo (ms)')
title(sprintf('Evoluzione temporale in x = %.1f', x(idx_mid)), Interpreter="latex")
axis([0, T, y_min y_max])
legend(Location="best")

for k = 1:length(t)
    addpoints(h1, t(k), V(idx_mid, k))
    addpoints(h2, t(k), W(idx_mid, k))
    pause(0.02)
end


figure; % plot evoluzione temporale 

plotV = plot(x, V(:, 1), 'b', LineWidth=2, DisplayName='potenziale');
hold on 
plotW = plot(x, W(:, 1), 'r', LineWidth=2, DisplayName='recovery');

ylim([y_min y_max])
legend(Location="best")
grid on
hold on 
xlabel('Intervallo spaziale', Interpreter='latex')

for j = 1:m+1
    set(plotV, 'YData', V(:, j))
    set(plotW, 'YData', W(:, j)) 
    title(sprintf('t = %.2f ms', t(j)), Interpreter="latex");
    drawnow
end




%% Simulazione nel caso in cui abbia due stimoli nei due estremi della fibra 
clear

sigma = 1e-3;
b = 5;
c = 1;
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

v0 = @(x) zeros(length(x), 1);
w0 = @(x) zeros(length(x), 1);

m = 800; % tempo
n = 1000; % spazio

T = 50; 

Dt = T/m;
h = 1/(n+1);

t = linspace(0, T, m+1)';
x = linspace(0, 1, n+2)';

V = zeros(n+2, m+1); % potenzale 
V(:, 1) = v0(x); 

W = zeros(n+2, m+1); % recovery
W(:, 1) = w0(x);

A = -(1/h^2) * spdiags(repmat([1, -2, 1], n+2, 1), -1:1, n+2, n+2);
A(1, 2) = -2/(h^2); A(n+2, n+1) = -2/(h^2); % valori di correzione per matrice A
A = (1/Dt).*speye(n+2) + sigma*A; 
dA = decomposition(A); % per risolvere il sistema lineare in modo più efficiente

% UNICA VARIAZIONE RISPETTO A PRIMA 
I_app = @(t, x) 2*((0 <= t & t <= 1) & ( (0 <= x & x <= 0.04) | (0.96 <= x & x <= 1 ) )); 

for j=2:m+1
    F = V(:, j-1)/Dt + b * V(:, j-1) .* ( V(:, j-1)-beta ) .* ( delta - V(:, j-1) ) - c*W(:, j-1) + I_app(t(j), x);

    V(:, j) = dA\F;
    W(:, j) = (Dt / (1+Dt*e*gamma)) * (W(:, j-1)./Dt + e*V(:, j-1));
end

figure;
imagesc(t, x, V); axis xy
xlabel('tempo (ms)')
ylabel('spazio')
title('Potenziale', Interpreter='latex')

figure;
mesh(t, x, V);
xlabel('tempo (ms)')
ylabel('spazio')
title('Potenziale', Interpreter='latex')


figure;
imagesc(t, x, W); axis xy
xlabel('tempo (ms)')
ylabel('spazio')
title('Recovery', Interpreter='latex')


figure;
mesh(t, x, W); axis xy
xlabel('tempo (ms)')
ylabel('spazio')
title('Recovery', Interpreter='latex')

figure; % plot evoluzione temporale nel punto x = 1/4

y_min = min(min(V(:)), min(W(:))) - 0.05;
y_max = max(max(V(:)), max(W(:))) + 0.05;
idx_x = ceil((n+2)/4); % indice punto x

h1 = animatedline('LineWidth', 1.5, Color='blue', DisplayName='potenziale');
h2 = animatedline('LineWidth', 1.5, Color='red', DisplayName='recovery');

grid on
hold on
xlabel('tempo (ms)')
title(sprintf('Evoluzione temporale in x = %.2f', x(idx_x)), Interpreter="latex")
axis([0, T, y_min y_max])
legend(Location="best")

for k = 1:length(t)
    addpoints(h1, t(k), V(idx_x, k))
    addpoints(h2, t(k), W(idx_x, k))
    pause(0.02)
end


figure; % plot evoluzione temporale 

plotV = plot(x, V(:, 1), 'b', LineWidth=2, DisplayName='potenziale');
hold on 
plotW = plot(x, W(:, 1), 'r', LineWidth=2, DisplayName='recovery');

ylim([y_min y_max])
legend(Location="best")
grid on
hold on 
xlabel('Intervallo spaziale', Interpreter='latex')

for j = 1:m+1
    set(plotV, 'YData', V(:, j))
    set(plotW, 'YData', W(:, j)) 
    title(sprintf('t = %.2f ms', t(j)), Interpreter="latex");
    drawnow
end






