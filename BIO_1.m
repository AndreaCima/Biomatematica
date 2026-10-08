clear; clc; close all; 

%% punto 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulare il sistema di reazione enzimatica con i parametri descritti nel
% modello del pd biomat_1.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametri del modello
k1 = 4e6;
km1 = 25;
k2 = 15; 
e0 = 1e-8;
tspan = [0, 100];
s0 = ( km1 + k2 ) / k1; 
Km = ( km1 + k2 ) / k1;

% definizione sistema di equazioni differenziali
y0 = [s0; 0];
options=odeset('RelTol',5.e-13 ,'AbsTol',1.e-13,'InitialStep',1.e-5,'MaxStep',5);

odefun = @(t, y) [km1*y(2) - k1*y(1)*(e0-y(2)); k1*y(1)*(e0-y(2)) - (km1+k2)*y(2)]; % s := y(1);  c := y(2)
[t, y] = ode15s(odefun, tspan, y0, options);

% normalizzazione
s = y(:, 1) ./ s0; 
c = y(:, 2) ./ e0;
% Le soluzioni riportate sopra sono quelle che poi saranno considerate come
% quelle "esatte", rispetto alle quali valuteremo errori di approssimazione
% commessi con altre approssimazioni

figure; 
plot(t, s, '-o', LineWidth=2, DisplayName='Substrato')
hold on
grid on
plot(t, c, '-d', LineWidth=2, DisplayName= 'Composto')
xlabel('tempo(s)')
legend(Location="best")
title('Evoluzione in tempo delle concentrazioni di substrato e composto', Interpreter='latex')

figure; 
plot(s, c, '-o', LineWidth=2)
title('Spazio delle fasi', Interpreter='latex')
xlabel('Substrato')
ylabel('Composto')
grid on

%% punto 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analizzando i plot dell'evoluzione di composto e substrato individuo in
% modo indicativo le due regioni di regimi transitorio e quasi stazionario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% regime transitorio
figure;
div = 36; % separazione tra le due regioni, individuata in modo indicativo
t_tra = t(1:div);
c_tra = c(1:div);

plot(t_tra, c_tra, '-o', LineWidth=2)
ylim([0, 0.5])
grid on
title('Composto: fase transitoria', Interpreter='latex')
xlabel('tempo(s)')

% regime quasi stazionario
figure; 
t_qs = t(div+1:end);
c_qs = c(div+1:end);
plot(t_qs, c_qs, '-o', LineWidth=2)
grid on
ylim([0, 0.5])
title('Composto: fase quasi-stazionaria', Interpreter='latex')
xlabel('tempo(s)')

%% punto 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulazione all'equilibrio e quasi stazionaria, confronto con soluzione
% completa. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uso l'approssimazione all'equilibrio c = e0 * s / (s+KM) e la sostituisco
% nell'ODE iniziale, dunque rimane solo l'incognita s

KM = km1/k1;
odefun_eq = @(t, s) km1 * (e0*s ./ (s + KM)) - k1*s*(e0 - (e0*s ./ (s + KM)));

%%%%%%%%%%%%%%%%%%% APPROSSIMAZIONE ALL'EQUILIBRIO %%%%%%%%%%%%%%%%%%%%%%%%
% in ode15s uso come tspan il vettore t che avevo già trovato nel calcolare
% la soluzione "esatta" nel punto 1
[t_approx_eq, s_approx_eq] = ode15s(odefun_eq, t, s0, options);

% devo normalizzare dividendo per e0, quindi non lo moltiplico davanti
c_approx_eq = s_approx_eq ./ (s_approx_eq + KM);

%%%%%%%%%%%%%%%%%%%%% APPROSSIMZIONE QUASI STAZIONARIA %%%%%%%%%%%%%%%%%%%%
% passaggi analoghi a quelli per l'approssimazione all'equilibrio
% In questo caso sostituisco nell'equazione iniziale c = e0 * s / (s+Km),
% dunque rimane solo l'incognita s 

odefun_qs =@(t, s) km1 * (e0*s ./ (s + Km)) - k1*s*(e0 - (e0*s ./ (s + Km)));

[t_approx_qs, s_approx_qs] = ode15s(odefun_qs, t, s0, options);

% devo normalizzare dividendo per e0, quindi non lo moltiplico davanti
c_approx_qs = s_approx_qs ./ (s_approx_qs + Km);

figure; 
plot(t, c,  LineWidth=2, DisplayName='Soluzione esatta')
hold on 
plot(t_approx_eq, c_approx_eq, '-o',  LineWidth=1, DisplayName="Approssimazione all'equilibrio")
plot(t_approx_qs, c_approx_qs, '-d',  LineWidth=1, DisplayName='Approssimazione quasi-stazionaria')

xlabel('tempo(s)')

legend(Location= "best")
title('Composto', Interpreter='latex')
grid on

% Mi sembra sensato che nell'approssimazione all'equilibrio c_approx_eq sia
% costante dal momento che in questa approssimazione sto supponedo ds/dt=0.
% Quindi s è costante e dunque anche c lo è, dal momento che è definito in
% termini di s.


%% punto 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Approssimazione uniforme delle concentrazioni s, c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

odefun_unif = @(t, s) -k2*e0*s ./ (Km + s);
[t_unif, s_unif] = ode15s(odefun_unif, t, s0, options);
c_unif = e0*s_unif ./ (Km + s_unif) - e0*s0 ./ (Km + s0) .* exp(-(Km+s0).*k1.*t_unif);

% normalizzazione
s_unif = s_unif ./ s0;
c_unif = c_unif ./ e0;

figure; 
plot(t, s,  DisplayName='Substrato')
hold on
plot(t_unif, s_unif, 'o', DisplayName='Substrato: approssimazione')
plot(t, c, DisplayName='Composto')
plot(t_unif, c_unif, 'o', DisplayName='Composto:approssimazione')

grid on
legend()
title('Approssimazione uniforme', Interpreter='latex')
xlabel('tempo(s)')

%% punto 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Studiare errori relativi di approssimazione per substrato e composto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err_s = abs(s - s_unif) ./ norm(s, 2);
err_c = abs(c - c_unif) ./ norm(c, 2);

figure; 
semilogy(t, err_s, LineWidth=2, DisplayName='Substrato')
hold on
semilogy(t, err_c, LineWidth=2, DisplayName='Composto')
legend(Location="best")
grid on
title('Errore relativo di approssimazione uniforme', Interpreter='latex')

%% punto 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analizzo errori di approssimazione separatamente su regione transitoria e
% quasi stazionaria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure; 
semilogy(t(1:div), err_s(1:div), LineWidth=2, DisplayName='Substrato')
hold on
semilogy(t(1:div), err_c(1:div), LineWidth=2, DisplayName='Composto')
legend(Location="best")
title('Errore relativo: fase transitoria', Interpreter='latex')
grid on

figure; 
semilogy(t(div+1:end), err_s(div+1:end), LineWidth=2, DisplayName='Substrato')
hold on
semilogy(t(div+1:end), err_c(div+1:end), LineWidth=2, DisplayName='Composto')
legend(Location="best")
title('Errore relativo: fase quasi stazionaria', Interpreter='latex')
grid on


%% punto 7
figure; 
plot(s, c, LineWidth=2, DisplayName='Curva esatta')
hold on 
plot(s_unif, c_unif, '-d', DisplayName='Approssimazione uniforme')
legend(Location="best")
grid on
title('Spazio delle fasi', Interpreter='latex')
xlabel('Substrato')
ylabel('Composto')


