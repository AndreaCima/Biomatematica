function V = Nerst(T, z, ce, ci)

R = 8.31; % costante gas perfetti
F = 96.49; % costante Faraday

V = R*T/(z*F) * log(ce/ci);
end
