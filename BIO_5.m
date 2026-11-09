clear; clc; close all

Cm = 1; 
GNa = 120;
GK = 36;
GL = 0.3;
VNa = 115;
VK = -12;
VL = 10.6;

alpha_m = @(V) 0.1*(25-V).*(exp((25-V)/10)-1).^(-1);
alpha_h = @(V) 0.07*exp(-V/20);
alpha_n = @(V) 0.01*(10-V) .* (exp( (10-V)/10 ) -1 ).^(-1);

beta_m = @(V) 4*exp(-V/18);
beta_h = @(V) (exp( (30-V)/10 ) + 1).^(-1);
beta_n = @(V) 0.125*exp(-V/80);

V0 = 6.6;
m0 = 5.2934e-2;
h0 = 5.9611e-1;
n0 = 3.1768e-1;
Iapp = 0;

T = 20;

tspan = [0, T];

y0 = [V0 m0 h0 n0]';

f = @(t, y) [(1/Cm) * (-(GNa*( y(2).^3 ) .* y(3) .* (y(1)-VNa) + GK*( y(4).^4 ).*( y(1)-VK ) + GL*( y(1)-VL ) ) + Iapp); ...
    alpha_m(y(1)).*(1-y(2)) - beta_m(y(1)).*y(2); ...
    alpha_h(y(1)).*(1-y(3)) - beta_h(y(1)).*y(3); ...
    alpha_n(y(1)).*(1-y(4)) - beta_n(y(1)).*y(4)];

[t, ya] = ode15s(f, tspan, y0);

V = ya(:, 1);
m = ya(:, 2);
h = ya(:, 3);
n = ya(:, 4);


tiledlayout(2, 2)
nexttile
plot(t, V, LineWidth=2)
grid on
xlabel('tempo (ms) ', Interpreter='latex')
title('Potenziale', Interpreter='latex')

nexttile
plot(t, m, LineWidth=2, DisplayName='m')
grid on
hold on
plot(t, h, LineWidth=2, DisplayName='h')
plot(t, n, LineWidth=2, DisplayName='n')
legend(Location="best")
xlabel('tempo (ms)')

nexttile
gNa = GNa * (m.^3) .* h;
gK = GK * n.^4;

plot(t, gNa, LineWidth=2, DisplayName='gNa')
grid on
hold on
plot(t, gK, LineWidth=2, DisplayName='gK')
title('Conduttanze', Interpreter='latex')
xlabel('tempo (ms)')
legend(Location="best")

nexttile 
INa = gNa .* (V-VNa);
IK = gK .*(V-VK);
IL = GL*(V-VL);

Iion = INa + IK + IL;


plot(t, INa, LineWidth=2, DisplayName='INa')
grid on 
hold on
plot(t, IK, LineWidth=2, DisplayName='IK')
plot(t, IL, LineWidth=2, DisplayName='IL')
plot(t, Iion, LineWidth=2, DisplayName='Iion')
legend(Location="best")



figure;
tau_m = 1./(alpha_m(V) + beta_m(V));
tau_h = 1./(alpha_h(V) + beta_h(V));
tau_n = 1./(alpha_n(V) + beta_n(V));

plot(V, tau_m, LineWidth=2, DisplayName='$\tau_m$')
grid on 
hold on
plot(V, tau_h, LineWidth=2, DisplayName='$\tau_h$')
plot(V, tau_n, LineWidth=2, DisplayName='$\tau_n$')
legend(Location="best", Interpreter="latex")
title('Costanti di tempo', Interpreter='latex')

figure;

m_inf = alpha_m(V)./(alpha_m(V) + beta_m(V));
h_inf = alpha_h(V)./(alpha_h(V) + beta_h(V));
n_inf = alpha_n(V)./(alpha_n(V) + beta_n(V));

plot(V, m_inf, LineWidth=2, DisplayName='$m_\infty$')
grid on
hold on
plot(V, h_inf, LineWidth=2, DisplayName='$h_\infty$')
plot(V, n_inf, LineWidth=2, DisplayName='$n_\infty$')
legend(Location="best", Interpreter="latex")
title('Valori qusi stazionari', Interpreter='latex')


