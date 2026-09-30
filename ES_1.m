clear; clc; close all; 

% parametri del modello
k1 = 4e6;
km1 = 25;
k2 = 15; 
e0 = 1e-8;
tspan = [0, 100];

% definizione sistema di equazioni differenziali
s0 = ( km1 + k2 ) / k1; 
y0 = [s0; 0];
options=odeset('RelTol',5.e-13 ,'AbsTol',[1.e-13 1.e-13],'InitialStep',1.e-5,'MaxStep',5);

odefun = @(t, y) [km1*y(2) - k1*y(1)*(e0-y(2)); k1*y(1)*(e0-y(2)) - (km1+k2)*y(2)]; % s := y(1);  c := y(2)
[t, y] = ode15s(odefun, tspan, y0, options);

% normalizzazione
s = y(:, 1) ./ s0; 
c = y(:, 2) ./ e0;

figure; 
plot(t, s, '-o', LineWidth=2, DisplayName='substrato')
hold on
grid on
plot(t, c, '-d', LineWidth=2, DisplayName= 'composto')
xlabel('tempo(s)')
legend(Location="best")

figure; 
plot(s, c, '-o', LineWidth=2)
title('Spazio delle fasi')
xlabel('substrato')
ylabel('composto')
grid on

% regime transitorio
figure;
div = 36;
t_tra = t(1:div);
c_tra = c(1:div);

plot(t_tra, c_tra, '-o', LineWidth=2)
ylim([0, 0.5])
grid on
title('composto: fase transitoria')
xlabel('tempo(s)')


figure; 
t_qs = t(div+1:end);
c_qs = c(div+1:end);
plot(t_qs, c_qs, '-o', LineWidth=2)
grid on
ylim([0, 0.5])
title('composto: fase quasi-stazionaria')
xlabel('tempo(s)')
legend()

% approx all'equilibrio
s_qs_real = y(div+1:end, 1);
KM = km1/k1;
c_app_eq = s_qs_real ./ (s_qs_real + KM);

% approssimazione quasi stazionaria
Km = (km1 + k2) / k1;
c_app_qs = s_qs_real ./ (s_qs_real + Km);

figure; 
plot(t_qs, c_qs,  LineWidth=2, DisplayName='fase quasi stazionaria')
hold on 
plot(t_qs, c_app_eq, '--',  LineWidth=1, DisplayName="approssimazione all'equilibrio")
plot(t_qs, c_app_qs, '--',  LineWidth=1, DisplayName='approssimazione quasi-stazionaria')

xlabel('tempo(s)')
ylabel('composto')

legend(Location= "best")
grid on

