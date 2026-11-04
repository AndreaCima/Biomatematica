clear; clc; close all; 
 
b = 5; 
c = 1; 
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

T = 100;
I_app = 0.0;

% v0 = 0; w0 = 0;
% v0 = 0.1; w0 = 0;
v0 = 0.6; w0 = 0;

y0 = [v0; w0];
tspan = [0 T];

f = @(t, y) [b*y(1).*(y(1)-beta).*(delta-y(1)) - c*y(2) + I_app; ...
    e*(y(1)-gamma*y(2))];

[t, ye] = ode15s(f, tspan, y0);

v = ye(:, 1);
w = ye(:, 2);

subplot(1, 2, 1)
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');

axis([0, T, min(ye(:)), max(ye(:))])
grid on
title('Evoluzione temporale', Interpreter='latex')
xlabel('tempo (s)', Interpreter='latex')
legend('potenziale', 'recovery', Location='best')

subplot(1, 2, 2)
h3 = animatedline('LineWidth', 1.5);
hold on
g = @(x) x./gamma;
f = @(x) (b/c)*x.*( x-beta ) .* ( delta-x ) + I_app;
% plot nullclines
fplot(f, [min(ye(:)), max(ye(:))])
fplot(g, [min(ye(:)), max(ye(:))])

F =@(x, y) b*x.*( x-beta ) .* ( delta-x ) -c*y+I_app;
G=@(x, y) e*(x-gamma*y);

[X, Y] = meshgrid(linspace(min(ye(:)), max(ye(:)), 30));
U = F(X, Y);
V = G(X, Y);
quiver(X, Y, U, V, Color='blue')

axis([min(ye(:)), max(ye(:)), min(ye(:)), max(ye(:))])
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
% considero un impulso di corrente applicata nell'intervallo [0, T_stim]

b = 5; 
c = 1; 
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;
alpha = 1; 

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

[t, ye] = ode15s(f, tspan, y0);

v = ye(:, 1);
w = ye(:, 2);

figure; 

subplot(1, 2, 1)
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');

axis([0, T, min(ye(:)), max(ye(:))])
grid on
title('Evoluzione temporale', Interpreter='latex')
xlabel('tempo (s)', Interpreter='latex')
legend('potenziale', 'recovery', Location='best')

subplot(1, 2, 2)
h3 = animatedline('LineWidth', 1.5);


axis([min(ye(:)), max(ye(:)), min(ye(:)), max(ye(:))])
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
