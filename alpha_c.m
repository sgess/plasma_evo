function alpha = alpha_c(n,T_e)
% n in cm^-3, T in eV
% alpha in cm^3/s

alpha = 8.75e-27*(T_e).^(-9/2).*n;