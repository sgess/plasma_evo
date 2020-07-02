function plasma_sim = init_sim_iTherm(plasma_sim)

SI_consts;
plasma_sim.electron_mass = SI_em; % kg
plasma_sim.ion_mass = plasma_sim.ion_mass_ratio*SI_pm; % kg

plasma_sim.dr = plasma_sim.r_ax(2) - plasma_sim.r_ax(1); % cm
plasma_sim.dt = plasma_sim.t_ax(2) - plasma_sim.t_ax(1); % cm

plasma_sim.nu_ii_prefactor = (SI_em/plasma_sim.ion_mass)^(1/2)/sqrt(2); % no unit
plasma_sim.nu_ei_prefactor = 2*(SI_em/plasma_sim.ion_mass); % no unit
plasma_sim.k_prefactor = 3.2*(4*pi*SI_eps0)^2*3/(4*sqrt(2*pi*SI_em)*SI_e^(3/2))/100;
plasma_sim.ki_prefactor = 3.9*(4*pi*SI_eps0)^2*3/(4*sqrt(2*pi*plasma_sim.ion_mass)*SI_e^(3/2))/100;
plasma_sim.D_ii_prefactor = 100^2*SI_e/plasma_sim.ion_mass; % no unit

plasma_sim.density = zeros(plasma_sim.n_r,plasma_sim.n_t); % plasma density, cm^-3
plasma_sim.neutral = zeros(plasma_sim.n_r,plasma_sim.n_t); % neutral density, cm^-3
plasma_sim.T_ions = zeros(plasma_sim.n_r,plasma_sim.n_t); % ion temperature, eV
plasma_sim.T_eles = zeros(plasma_sim.n_r,plasma_sim.n_t); % electron temperature, eV
plasma_sim.T_neut = zeros(plasma_sim.n_r,plasma_sim.n_t); % neutral temperature, eV
plasma_sim.Nu_ees = zeros(plasma_sim.n_r,plasma_sim.n_t); % electron-electron collision freuqency, s^-1
plasma_sim.Nu_i0s = zeros(plasma_sim.n_r,plasma_sim.n_t); % electron-ion collision frequency, s^-1
plasma_sim.D_ambi = zeros(plasma_sim.n_r,plasma_sim.n_t); % ambipolar diffusion, cm^2 s^-1
plasma_sim.alpha3 = zeros(plasma_sim.n_r,plasma_sim.n_t); % 3-body reco rate, cm^3 s^-1
plasma_sim.k_therm = zeros(plasma_sim.n_r,plasma_sim.n_t); % plasma thermal conductivity, cm^-1 s^-1
plasma_sim.ki_therm = zeros(plasma_sim.n_r,plasma_sim.n_t); % plasma thermal conductivity, cm^-1 s^-1

% Coulomb logarithm
lambda = lambda_ei(plasma_sim.n_init,plasma_sim.T_ele_init);

% Electron-Electron collision frequency s^-1
nu_ee = nu_physics(plasma_sim.n_init,plasma_sim.T_ele_init);

% Ion-Ion collision frequency s^-1
nu_ii = plasma_sim.nu_ii_prefactor*(plasma_sim.T_ele_init./plasma_sim.T_ion_init).^(3/2).*nu_ee;

% Ion thermal velocity cm/s
v_ion = thermal_velocity(plasma_sim.T_ion_init,plasma_sim.ion_mass);

% Ion-neutral cross section cm^2
sigma = RbCross(v_ion);

% Neutral density cm^-3
n_n = plasma_sim.n_0 - plasma_sim.n_init + plasma_sim.n_min;

% Ion-neutral collision frequency s^-1
nu_0 = nu_io(n_n,sigma,v_ion);

% Ion Diffusion rate cm^2 s^-1
D_ii = plasma_sim.D_ii_prefactor*plasma_sim.T_ion_init./(nu_ii+nu_0);

% Ambipolar Diffusion rate cm^2 s^-1
D_a = (1+plasma_sim.T_ele_init./plasma_sim.T_ion_init).*D_ii;

% Thermal Conductivity cm^-1 s^-1
k_therm = plasma_sim.k_prefactor.*plasma_sim.T_ele_init.^(5/2)./lambda;
ki_therm = plasma_sim.ki_prefactor.*plasma_sim.T_ion_init.^(5/2)./lambda;


% recombination coeff cm^3 s^-1
alpha3 = alpha_c(plasma_sim.n_init,plasma_sim.T_ele_init);

plasma_sim.density(:,1) = plasma_sim.n_init;
plasma_sim.neutral(:,1) = n_n;
plasma_sim.T_ions(:,1) = plasma_sim.T_ion_init;
plasma_sim.T_eles(:,1) = plasma_sim.T_ele_init;
plasma_sim.T_neut(:,1) = plasma_sim.T_neut_init;
plasma_sim.Nu_ees(:,1) = nu_ee;
plasma_sim.Nu_i0s(:,1) = nu_0;
plasma_sim.D_ambi(:,1) = D_a;
plasma_sim.alpha3(:,1) = alpha3;
plasma_sim.k_therm(:,1) = k_therm;
plasma_sim.ki_therm(:,1) = ki_therm;
