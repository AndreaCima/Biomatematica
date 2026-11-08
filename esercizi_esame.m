%% 2019 esercizio 3 
clear; clc; close all; 

% parametri
b = 5; 
beta = 0.1; 
epsilon = 0.1; 
gamma = 0.25; 

% diverse scelte di dati iniziali
% v0 = 0.2; w0 = 0.1; 
% v0 = 0.3; w0 = 0.1;
% v0 = -0.1; w0 = -0.1; 
v0 = 0.13; w0 = 0.1;

T = 300; 

y0 = [v0 w0]';
tspan = [0, T];

% notazione v = y(1); w = y(2)
f = @(t, y) [b * y(1) .* ( y(1)-beta ) .* ( 1-y(1) ) - y(1).*y(2); ...
    epsilon * ( y(1)-gamma*y(2) ) ];

[t, ya] = ode15s(f, tspan, y0);

min_ya = min(ya(:)) - 0.05; % aggiungo un piccolo margine per i plot
max_ya = max(ya(:)) + 0.05;

v = ya(:, 1);
w = ya(:, 2);

f_null = @(x) b.*(x-beta).*(1-x);

g_null = @(x) x/gamma; 

F = @(x, y) b.*x.*(x-beta).*(1-x) - x.*y;
G = @(x, y) epsilon*(x-gamma*y);

[X, Y] = meshgrid(linspace(min_ya, max_ya, 20));

% per plot campo vettoriale con quiver
U = F(X, Y);
V = G(X, Y);


tiledlayout(1, 2)

nexttile
h1 = animatedline(LineWidth=2, Color='b');
hold on
h2 = animatedline(LineWidth=2, Color='r');
title('Evoluzione temporale', Interpreter='latex')
grid on
xlabel('tempo (s)')
legend('potenziale', 'recovery', Location='best')
axis([0 T min_ya max_ya])

nexttile
h3 = animatedline(LineWidth=2, Color='k');
hold on
grid on
fplot(f_null, [min_ya max_ya], Color='b')
xline(0, Color = 'b') % f si annulla anche per v = 0
fplot(g_null, [min_ya max_ya], Color='r')

quiver(X, Y, U, V, Color='b')
title('Spazio delle fasi', Interpreter='latex')
axis([min_ya max_ya min_ya max_ya])

for k = 1:length(t)
    addpoints(h1, t(k), v(k))
    addpoints(h2, t(k), w(k))
    addpoints(h3, v(k), w(k))
    drawnow limitrate
    pause(0.1)
end
drawnow


