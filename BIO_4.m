clear; clc; close all; 

 
b = 5; 
c = 1; 
beta = 0.1;
delta = 1;
gamma = 0.25;
e = 0.1;

T = 100;
I_app = 0.5;
v0 = 0; w0 = 0;
% v0 = 0.1; w0 = 0;
% v0 = 0.6; w0 = 0;

y0 = [v0; w0];
tspan = [0 T];
f = @(t, y) [b*y(1).*(y(1)-beta).*(delta-y(1)) - c*y(2) + I_app; ...
    e*(y(1)-gamma*y(2))];

[t, ye] = ode15s(f, tspan, y0);
v = ye(:, 1);
w = ye(:, 2);

tiledlayout(1, 2)

title('FHN in 0D', Interpreter='latex')

nexttile
h1 = animatedline('LineWidth', 1.5, Color='blue');
h2 = animatedline('LineWidth', 1.5, Color='red');
axis([0, T, min(ye(:)), max(ye(:))])
grid on
title('Evoluzione temporale', Interpreter='latex')

nexttile
h3 = animatedline('LineWidth', 1.5);
axis([min(ye(:)), max(ye(:)), min(ye(:)), max(ye(:))])
grid on
title('Spazio delle fasi', Interpreter='latex')

for k = 1:length(t)
    addpoints(h1, t(k), v(k))
    addpoints(h2, t(k), w(k))
    drawnow

    addpoints(h3, v(k), w(k))
    drawnow

    pause(0.2)
end

