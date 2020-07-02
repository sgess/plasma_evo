function v = thermal_velocity(T,m)
% T in eV, m in kg
% v in cm/s

SI_consts;

%v = 100*SI_c * sqrt(3*T_i/(1e6*A_i*SI_pM));
v = 100 * sqrt(3*SI_e*T/m);