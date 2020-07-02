%% plasma_sim struct will contain all of the quantities needed for the calculation

plasma_sim = struct();

plasma_sim.name = 'example';

plasma_sim.r_pipe = 2; % cm
plasma_sim.n_r = 1001; % number of radial points
plasma_sim.r_ax = linspace(0,plasma_sim.r_pipe,plasma_sim.n_r)'; % cm

plasma_sim.t_max = 1E-4; % seconds
plasma_sim.n_t = 1001; % number of time points
plasma_sim.t_ax = linspace(0,plasma_sim.t_max,plasma_sim.n_t); % seconds

plasma_sim.n_0 = 1.81E14; % cm^-3
plasma_sim.n_min = 10; % cm^-3

plasma_sim.T_e0 = 0.41; % eV
plasma_sim.T_i0 = 0.041; % eV

plasma_sim.ion_mass_ratio = 84.85; % in units of proton mass
plasma_sim.ionization_energy = 4.18; % eV

%% Set up radial plasma density profile initial condition

% Use results of Gabor's sim as input for density profile
load('plasma_ICs_All.mat');
energy_string = '_135mJ';
order_string  = '_0';
plasma_profile = laser_input.(['plasma_profile' energy_string order_string]);
prof = interp1(laser_input.rr/10,plasma_profile,plasma_sim.r_ax,'linear','extrap');
plasma_sim.n_init = (plasma_sim.n_0-plasma_sim.n_min)*prof+plasma_sim.n_min;

% Or create a smooth function that approximate's profile from Gabor's sim

%r0 = 0.1555; % cm for 135 mJ, r_fwhm = 0.216 mm
%ord = 10; % cm for 135 mJ

%r0 = 0.1468; % cm for 95 mJ
%ord = 10; % cm for 95 mJ

%r0 = 0.1210; % cm for 40 mJ
%ord = 7; % cm for 40 mJ

%plasma_sim.n_init = (plasma_sim.n_0-plasma_sim.n_min)*exp(-((plasma_sim.r_ax).^2/(2*r0^2)).^ord) + plasma_sim.n_min;

%% Set up radial temperature profile initial condition

% Use results of Gabor's sim as input for temperature profile
energy_profile = laser_input.(['energy_profile' energy_string order_string]);
eng = interp1(laser_input.rr/10,energy_profile,plasma_sim.r_ax,'linear','extrap');
plasma_sim.T_ele_init = (2/3)*eng;
plasma_sim.T_ele_init(plasma_sim.T_ele_init < plasma_sim.T_i0) = plasma_sim.T_i0; 

% Or create a smooth function that approximate's profile from Gabor's sim
%plasma_sim.T_ele_init = (plasma_sim.T_e0-plasma_sim.T_i0)*exp(-((plasma_sim.r_ax).^2/(2*r0^2)).^ord)+plasma_sim.T_i0;

% Ion and neutral temperatures start spatially uniform at vapor temp
plasma_sim.T_ion_init = plasma_sim.T_i0*ones(size(plasma_sim.r_ax));
plasma_sim.T_neut_init = plasma_sim.T_i0*ones(size(plasma_sim.r_ax));

%% Plot initial conditions

figure();
plot(plasma_sim.r_ax,plasma_sim.n_init,'b','linewidth',2);
xlabel('R [cm]');
ylabel('Density [cm^{-3}]');
set(gca,'fontsize',16);

figure();
plot(plasma_sim.r_ax,plasma_sim.T_ele_init,'b','linewidth',2);
xlabel('R [cm]');
ylabel('Temperature [eV]');
set(gca,'fontsize',16);

%% Initialize simulation

plasma_sim = init_sim_iTherm(plasma_sim);

figure();
semilogy(plasma_sim.r_ax,plasma_sim.D_ambi(:,1),'b','linewidth',2);
xlabel('R [mm]');
ylabel('Diffusion constant (logscale) D_a [cm^2/s]');
set(gca,'fontsize',16);

figure();
semilogy(plasma_sim.r_ax,plasma_sim.Nu_ees(:,1),'b','linewidth',2);
xlabel('R [mm]');
ylabel('Electron-Electron collision frequency (logscale) \nu [s^{-1}]');
set(gca,'fontsize',16);

figure();
plot(plasma_sim.r_ax,plasma_sim.k_therm(:,1),'b','linewidth',2);
xlabel('R [mm]');
ylabel('\kappa_e [cm^{-1} s^{-1}]');
set(gca,'fontsize',16);

figure();
plot(plasma_sim.r_ax,plasma_sim.ki_therm(:,1),'b','linewidth',2);
xlabel('R [mm]');
ylabel('\kappa_i [cm^{-1} s^{-1}]');
set(gca,'fontsize',16);

figure();
plot(plasma_sim.r_ax,plasma_sim.alpha3(:,1),'b','linewidth',2);
xlabel('R [mm]');
ylabel('\alpha_3 [cm^3 s^{-1}]');
set(gca,'fontsize',16);

%% Run simulation

% if using a fixed boundary, the wall will be held at a specific temp 
% for both ions and electrons
plasma_sim.fixed_boundaries = true;
plasma_sim.Te_edge = 0.13;
plasma_sim.Ti_edge = plasma_sim.T_i0;


% if using sheath boundary, the B.C. for electrons is 
% Neumann with dT/dr = 0 at the boundary and 
% ions have the temperature of the boundary
plasma_sim.sheath_boundaries = false;
%plasma_sim.Te_edge = plasma_sim.T_i0;
%plasma_sim.Ti_edge = plasma_sim.T_i0;

% if using dynamic boundary, the wall will be heated by electrons and ions
% but the B.C. is still fixed boundary and the boundary temperature changes with time
plasma_sim.heated_boundaries = false;
%plasma_sim.Te_edge = plasma_sim.T_i0;
%plasma_sim.Ti_edge = plasma_sim.T_i0;
%plasma_sim.k_edge = 6.0E-11;


% run sim
plasma_sim = run_sim_bc_choice(plasma_sim);

%% Add experimental data to struct for comparison

awk_dat = load('data.txt');

plasma_sim.data.time = awk_dat(1,:)/1e9; % s
plasma_sim.data.high_eng = awk_dat(2,:); % cm^-3
plasma_sim.data.high_err = awk_dat(3,:); % cm^-3
plasma_sim.data.med_eng = awk_dat(4,:); % cm^-3
plasma_sim.data.med_err = awk_dat(5,:); % cm^-3
plasma_sim.data.low_eng = awk_dat(6,:); % cm^-3
plasma_sim.data.low_err = awk_dat(7,:); % cm^-3

% save simulation and data to single .mat file
% save([plasma_sim.name '.mat'],'plasma_sim');

%% Make some plots

figure();
mesh(1e6*plasma_sim.t_ax,plasma_sim.r_ax,plasma_sim.density/1e14); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
set(gca,'fontsize',14);

figure();
plot(plasma_sim.r_ax,plasma_sim.density(:,1)/1e14,'b',...
     plasma_sim.r_ax,plasma_sim.density(:,11)/1e14,'r',...
     plasma_sim.r_ax,plasma_sim.density(:,101)/1e14,'g',...
     plasma_sim.r_ax,plasma_sim.density(:,1001)/1e14,'k','linewidth',2);
xlabel('R [mm]');
ylabel('n [10^{14} cm^{-3}]');
legend({'t = 0 \mus','t = 1 \mus','t = 10 \mus','t = 100 \mus'});
set(gca,'fontsize',14);

% add a 1 ps offset to first time so that it is visible on log plot
model_t_plot = plasma_sim.t_ax; % s
model_t_plot(1) = 1e-12; % s
data_t_plot = plasma_sim.data.time;
data_t_plot(1) = 1e-12; % s

figure();
plot(1e6*model_t_plot,plasma_sim.density(1,:),'r','linewidth',2);
hold on;
errorbar(1e6*data_t_plot,plasma_sim.data.high_eng,plasma_sim.data.high_err,'bo','linewidth',2);
%errorbar(1e6*data_t_plot,plasma_sim.data.med_eng,plasma_sim.data.med_err,'bo','linewidth',2);
%errorbar(1e6*data_t_plot,plasma_sim.data.low_eng,plasma_sim.data.low_err,'bo','linewidth',2);

set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel('t [\mus]');
ylabel('n [10^{14} cm^{-3}]');
legend({'Model','Data'},'location','southwest');
set(gca,'fontsize',14);

figure();
plot(1e6*model_t_plot,plasma_sim.density(1,:),'r','linewidth',2);
hold on;
errorbar(1e6*data_t_plot,plasma_sim.data.high_eng,plasma_sim.data.high_err,'bo','linewidth',2);
%errorbar(1e6*data_t_plot,plasma_sim.data.med_eng,plasma_sim.data.med_err,'bo','linewidth',2);
%errorbar(1e6*data_t_plot,plasma_sim.data.low_eng,plasma_sim.data.low_err,'bo','linewidth',2);

set(gca,'xscale','log');
xlabel('t [\mus]');
ylabel('n [10^{14} cm^{-3}]');
legend({'Model','Data'},'location','southwest');
set(gca,'fontsize',14);

%%

figure();
pcolor(1e6*plasma_sim.t_ax,plasma_sim.r_ax,plasma_sim.density/1e14); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
set(gca,'fontsize',14);
colorbar;

figure();
pcolor(1e6*plasma_sim.t_ax,plasma_sim.r_ax,plasma_sim.T_eles); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
caxis([0 0.42]);
set(gca,'fontsize',14);
colorbar;

figure();
pcolor(1e6*plasma_sim.t_ax,plasma_sim.r_ax,plasma_sim.T_ions); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
caxis([0 0.42]);
set(gca,'fontsize',14);
colorbar;

figure();
pcolor(1e6*plasma_sim.t_ax,plasma_sim.r_ax,plasma_sim.T_neut); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
caxis([0 0.42]);
set(gca,'fontsize',14);
colorbar;

figure();
pcolor(1e6*plasma_sim.t_ax,plasma_sim.r_ax,log10(plasma_sim.density)); shading flat;
ylabel('R [mm]');
xlabel('t [\mus]');
zlabel('n [10^{14} cm^{-3}]');
set(gca,'fontsize',14);
colorbar;

figure();
plot(1e6*plasma_sim.t_ax,plasma_sim.T_eles(end,:),'r','linewidth',2);
xlabel('t [\mus]');
ylabel('Edge Temp. [eV]');
set(gca,'fontsize',14);