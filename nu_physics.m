function nu_ee = nu_physics(n,T)
% n in cm^-3, T in eV

% We use the definition from Callen 2.17
% pre_nu_ee_SI = 4*sqrt(2*pi)/(3*(4*pi)^2)*(SI_e^(5/2))/(SI_eps0^2*SI_em^(1/2))*1e6;
pre_nu_ee_SI = 2.9063e-06;

%lambda_ei = chen_coulomb_log(n,T);
l_ei = lambda_ei(n,T);

nu_ee = pre_nu_ee_SI.*n.*l_ei./T.^(3/2);