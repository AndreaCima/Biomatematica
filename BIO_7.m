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

gNa = GNa * (m.^3) .* h;
gK = GK * n.^4;

INa = gNa .* (V-VNa);
IK = gK .*(V-VK);
IL = GL*(V-VL);

Iion = INa + IK + IL;

min_V = min(V(:)) - 0.05; max_V = max(V(:)) + 0.05;
min_m = min(m(:)) - 0.05; max_m = max(m(:)) + 0.05;
min_h = min(h(:)) - 0.05; max_h = max(h(:)) + 0.05;
min_n = min(n(:)) - 0.05; max_n = max(n(:)) + 0.05;

min_gating = min([min_m, min_h, min_n]); 
max_gating = max([max_m, max_h, max_n]);

min_INa = min(INa(:)) - 0.05; max_INa = max(INa(:)) + 0.05;
min_IK = min(IK(:)) - 0.05; max_IK = max(IK(:)) + 0.05;
min_Iion = min(Iion(:)) - 0.05; max_Iion = max(Iion(:)) + 0.05;

min_I = min([min_INa, min_IK, min_Iion]);
max_I = min([max_INa, max_IK, max_Iion]);

figure; 
mesh(t, x, V)
xlabel('tempo')
ylabel('spazio')
title('Potenziale', Interpreter='latex')


% figure; 
% plotV = plot(x, V(:, 1), 'b', LineWidth=2, DisplayName='potenziale'); 
% hold on
% grid on 
% axis([0 1 min_V max_V])
% 
% for j = 1:n_t+1
%     set(plotV, 'YData', V(:, j))
%     title(sprintf('Potenziale, t = %.2f', t(j)), Interpreter="latex")
%     drawnow
%     % pause(0.05)
% end

figure;
tl = tiledlayout(1, 3);
nexttile 
plotV = plot(x, V(:, 1), 'b', LineWidth=2, DisplayName='potenziale'); 
grid on 
hold on
axis([0 1 min_V max_V])
title('Potenziale', Interpreter='latex')

nexttile
plotm = plot(x, m(:, 1), LineWidth=2, DisplayName='m');
hold on
grid on
ploth = plot(x, h(:, 1), LineWidth=2, DisplayName='h');
plotn = plot(x, n(:, 1), LineWidth=2, DisplayName='n');
legend(Location="best")
axis([0 1 min_gating max_gating])
title('Variabili di gating', Interpreter='latex')

nexttile
plot_INa = plot(x, INa(:, 1), LineWidth=2, DisplayName='INa');
hold on 
grid on
plot_IK = plot(x, IK(:, 1), LineWidth=2, DisplayName='IK');
plot_Iion = plot(x, Iion(:, 1), LineWidth=2, DisplayName='Iion');
legend(Location="best")
axis([0 1 min_I max_I])
title('Correnti', Interpreter='latex')

xlabel(tl, 'intervallo spaziale')
title(tl, sprintf('t = %.2f ms', t(1)), Interpreter='latex')
for k = 1:n_t+1
    % title(tiledlayout(1, 3), sprintf('t = %.2f', t(k)))
    set(plotV, 'YData', V(:, k))
    set(plotm, 'YData', m(:, k))
    set(ploth, 'YData', h(:, k))
    set(plotn, 'YData', n(:, k))
    set(plot_INa, 'YData', INa(:, k))
    set(plot_IK, 'YData', IK(:, k))
    set(plot_Iion, 'YData', Iion(:, k))
    tl.Title.String = sprintf('t = %.2f ms', t(k));
    drawnow
    pause(0.05)
end