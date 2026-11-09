clear; clc; close all; 
 
% parametri
b = 5; 
c = 1; 
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

T = 100;
I_app = 0.5;

% diverse scelte di condizioni inziali
v0 = 0; w0 = 0;
% v0 = 0.1; w0 = 0;
% v0 = 0.6; w0 = 0;

y0 = [v0; w0];
tspan = [0 T];

f = @(t, y) [b*y(1).*(y(1)-beta).*(delta-y(1)) - c*y(2) + I_app; ...
    e*(y(1)-gamma*y(2))];

[t, ya] = ode15s(f, tspan, y0);

v = ya(:, 1);
w = ya(:, 2);

min_ya = min(ya(:)) - 0.05; 
max_ya = max(ya(:)) + 0.05;

subplot(1, 2, 1)
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');

axis([0, T, min_ya, max_ya])
grid on
title('Evoluzione temporale', Interpreter='latex')
xlabel('tempo (s)', Interpreter='latex')
legend('potenziale', 'recovery', Location='best')

subplot(1, 2, 2)
h3 = animatedline('LineWidth', 1.5);
hold on

% plot nullclines
g_null = @(x) x./gamma;
f_null = @(x) (b*x.*( x-beta ) .* ( delta-x ) + I_app)/c;

fplot(f_null, [min_ya, max_ya])
fplot(g_null, [min_ya, max_ya])

F = @(x, y) b*x.*( x-beta ) .* ( delta-x ) -c*y + I_app;
G = @(x, y) e*(x-gamma*y);

[X, Y] = meshgrid(linspace(min_ya, max_ya, 30));
U = F(X, Y);
V = G(X, Y);
quiver(X, Y, U, V, Color='blue') % plot campo vettoriale

axis([min_ya, max_ya, min_ya, max_ya])
grid on
title('Spazio delle fasi', Interpreter='latex')
xlabel('v'); ylabel('w')

for k = 1:length(t)
    addpoints(h1, t(k), v(k))
    addpoints(h2, t(k), w(k))
    addpoints(h3, v(k), w(k))
    drawnow limitrate
    pause(0.1)
end
drawnow % per essere sicuro che l'animazione arrivi fino alla fine


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% considero un impulso di corrente applicata nell'intervallo [0, T_stim]

b = 5; 
c = 1; 
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;
alpha = 0.5; 

T = 100;
T_stim = 1;
I_app = @(t) alpha * (t <= T_stim);

v0 = 0; w0 = 0;
% v0 = 0.1; w0 = 0;
% v0 = 0.6; w0 = 0;

y0 = [v0; w0];
tspan = [0 T];

f = @(t, y) [b*y(1).*(y(1)-beta).*(delta-y(1)) - c*y(2) + I_app(t); ...
    e*(y(1)-gamma*y(2))];

[t, ya] = ode15s(f, tspan, y0);

min_ya = min(ya(:)) - 0.05;
max_ya = max(ya(:)) + 0.05;

v = ya(:, 1);
w = ya(:, 2);

figure; 
tiledlayout(1, 2)

nexttile
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');

axis([0, T, min_ya, max_ya])
grid on
title('Evoluzione temporale', Interpreter='latex')
xlabel('tempo (s)', Interpreter='latex')
legend('potenziale', 'recovery', Location='best')

nexttile
h3 = animatedline('LineWidth', 1.5);
hold on

% plot nullclines
g_null = @(x) x./gamma;
f_null = @(x) (b*x.*( x-beta ) .* ( delta-x ))/c;
fplot(f_null, [min_ya, max_ya])
fplot(g_null, [min_ya, max_ya])

% guardo il campo vettoriale F quando non ho più impulso, dal momento che
% questo campo vettoriale descrive il sistema per la maggior parte della
% simulazione
F = @(x, y) b*x.*( x-beta ) .* ( delta-x ) - c*y;
G = @(x, y) e*(x-gamma*y);

[X, Y] = meshgrid(linspace(min_ya, max_ya, 30));
U = F(X, Y);
V = G(X, Y);
quiver(X, Y, U, V, Color='blue') % plot campo vettoriale

axis([min_ya, max_ya, min_ya, max_ya])
grid on
title('Spazio delle fasi', Interpreter='latex')
xlabel('v'); ylabel('w')

for k = 1:length(t)
    addpoints(h1, t(k), v(k))
    addpoints(h2, t(k), w(k))
    addpoints(h3, v(k), w(k))
    drawnow limitrate
    pause(0.1)
end
drawnow % per essere sicuro che l'animazione arrivi fino alla fine
