function plasma_sim = run_sim_bc_choice(plasma_sim)

SI_consts;

n_r = plasma_sim.n_r;
n_t = plasma_sim.n_t;
dt = plasma_sim.dt;
dr = plasma_sim.dr;
q = dt/dr^2;

nu_ei_pre = plasma_sim.nu_ei_prefactor;
nu_ii_pre = plasma_sim.nu_ii_prefactor;
k_pre = plasma_sim.k_prefactor;
ki_pre = plasma_sim.ki_prefactor;
D_pre = plasma_sim.D_ii_prefactor;

% Generate plasma diffusion PDE matrix
M_d = matrix_D_neumann(n_r,plasma_sim.D_ambi(:,1),q);

% Generate electron temperature diffusion PDE matrix, for fixed walls
if plasma_sim.fixed_boundaries
    [M_t, bc_t] = matrix_T_input(n_r,plasma_sim.density(:,1),plasma_sim.k_therm(:,1),q,plasma_sim.Te_edge);
end

% Generate electron temperature diffusion PDE matrix, for sheath boundary (neumann bc)
if plasma_sim.sheath_boundaries
    M_t = matrix_T_neumann(n_r,plasma_sim.density(:,1),plasma_sim.k_therm(:,1),q);
    bc_t = 0;
end

% Generate electron temperature diffusion PDE matrix, for heated boundary
if plasma_sim.heated_boundaries
    denom = plasma_sim.k_edge/SI_boltz_J;
    ke_half = (1/2)*(plasma_sim.k_therm(end-1,1)+plasma_sim.k_therm(end,1));
    e_flux = ke_half*(plasma_sim.T_eles(end-1,1)-plasma_sim.T_eles(end,1))/dr;
    ki_half = (1/2)*(plasma_sim.ki_therm(end-1,1)+plasma_sim.ki_therm(end,1));
    i_flux = ki_half*(plasma_sim.T_ions(end-1,1)-plasma_sim.T_ions(end,1))/dr;
    plasma_sim.Te_edge = plasma_sim.Te_edge + (e_flux+i_flux)*dt/denom;
    plasma_sim.Ti_edge = plasma_sim.Te_edge;
    
    [M_t, bc_t] = matrix_T_input(n_r,plasma_sim.density(:,1),plasma_sim.k_therm(:,1),q,plasma_sim.Te_edge);
end

% Generate temperature diffusion PDE matrix
[M_i,bc_i] = matrix_T_input(n_r,plasma_sim.density(:,1),plasma_sim.ki_therm(:,1),q,plasma_sim.Ti_edge);


tic;
% Propagate forward
for i = 2:n_t
    
    
    % First, propagate plasma diffusion equation
    
    % Matrix in only applied to interior points
    N_now = plasma_sim.density(2:(end-1),i-1);
    
    % Density to be propagated forward. Recombination and B.C.s included at this step
    N_prop = N_now - dt*plasma_sim.alpha3(2:(end-1),i-1).*N_now.^2;
    
    % Forward propagation
    N_next = M_d\N_prop;
    
    % Update density and specify values at boundary
    plasma_sim.density(2:(end-1),i) = N_next;
    plasma_sim.density(1,i) = N_next(1);
    plasma_sim.density(end,i) = N_next(end);
    
    % Update neutral density
    plasma_sim.neutral(:,i) = plasma_sim.n_0 - plasma_sim.density(:,i) + plasma_sim.n_min;
    neutral_change = plasma_sim.neutral(:,i) - plasma_sim.neutral(:,i-1);
    
    % Update neutral temperature based on temperature of ions which have recombined.
    plasma_sim.T_neut(:,i-1) = (plasma_sim.T_neut(:,i-1).*plasma_sim.neutral(:,i-1) + plasma_sim.T_ions(:,i-1).*neutral_change)./plasma_sim.neutral(:,i);
    
    
    % Next, propagate thermal electron diffusion equation
    
    % Matrix in only applied to interior points
    T_now = plasma_sim.T_eles(2:(end-1),i-1);
    
    % Temperature to be propagated forward. Thermalization and B.C.s included at this step
    T_prop = T_now + dt*nu_ei_pre*plasma_sim.Nu_ees(2:(end-1),i-1).*(plasma_sim.T_ions(2:(end-1),i-1) - T_now) + bc_t;
    
    % Forward propagation
    T_next = M_t\T_prop;
    
    % Update electron temperature and specify values at boundary
    plasma_sim.T_eles(2:(end-1),i) = T_next;
    plasma_sim.T_eles(1,i) = T_next(1);
    
    if plasma_sim.sheath_boundaries
        plasma_sim.T_eles(end,i) = T_next(end);
    else
        plasma_sim.T_eles(end,i) = plasma_sim.Te_edge;
    end
    
    
    % Next, propagate thermal ion diffusion equation
    
    % Matrix in only applied to interior points
    Ti_now = plasma_sim.T_ions(2:(end-1),i-1);
    
    Ti_prop = Ti_now + dt*nu_ei_pre*plasma_sim.Nu_ees(2:(end-1),i-1).*(plasma_sim.T_eles(2:(end-1),i-1) - Ti_now) + ...
        2*dt*plasma_sim.Nu_i0s(2:(end-1),i-1).*(plasma_sim.T_neut(2:(end-1),i-1) - Ti_now) + bc_i;
    
    % Forward propagation
    Ti_next = M_i\Ti_prop;
    
    % Update ion temperature and specify values at boundary
    plasma_sim.T_ions(2:(end-1),i) = Ti_next;
    plasma_sim.T_ions(1,i) = Ti_next(1);
    plasma_sim.T_ions(end,i) = plasma_sim.Ti_edge;
    
    
    % Update neutral temperature
    plasma_sim.T_neut(:,i) = plasma_sim.T_neut(:,i-1) + 2*dt*plasma_sim.Nu_i0s(:,i-1).*(plasma_sim.T_ions(:,i-1) - plasma_sim.T_neut(:,i-1));
    
    
    
    % Update coefficients
    
    % Electron-Electron collision frequency s^-1
    plasma_sim.Nu_ees(:,i) = nu_physics(plasma_sim.density(:,i),plasma_sim.T_eles(:,i));
    
    % Ion-Ion collision frequency s^-1
    nu_ii = nu_ii_pre*(plasma_sim.T_eles(:,i)./plasma_sim.T_ions(:,i)).^(3/2).*plasma_sim.Nu_ees(:,i);
    
    % Ion thermal velocity cm/s
    v_ion = thermal_velocity(plasma_sim.T_ions(:,i),plasma_sim.ion_mass);
    
    % Ion-neutral cross section cm^2
    sigma = RbCross(v_ion);
    
    % Ion-neutral collision frequency s^-1
    plasma_sim.Nu_i0s(:,i) = nu_io(plasma_sim.neutral(:,i),sigma,v_ion);
    
    % Ion Diffusion rate cm^2 s^-1
    D_ii = D_pre*plasma_sim.T_ions(:,i)./(nu_ii+plasma_sim.Nu_i0s(:,i));
    
    % Ambipolar Diffusion rate cm^2 s^-1
    plasma_sim.D_ambi(:,i) = (1+plasma_sim.T_eles(:,i)./plasma_sim.T_ions(:,i)).*D_ii;
    
    % Coulomb logarithm
    lambda = lambda_ei(plasma_sim.density(:,i),plasma_sim.T_eles(:,i));
    
    % Thermal Conductivity cm^-1 s^-1
    plasma_sim.k_therm(:,i) = k_pre.*plasma_sim.T_eles(:,i).^(5/2)./lambda;
    plasma_sim.ki_therm(:,i) = ki_pre.*plasma_sim.T_ions(:,i).^(5/2)./lambda;
    
    % recombination coeff cm^3 s^-1
    plasma_sim.alpha3(:,i) = alpha_c(plasma_sim.density(:,i),plasma_sim.T_eles(:,i));
    

    
    % Finally, recalc matrices
    
    % Generate plasma diffusion PDE matrix
    M_d = matrix_D_neumann(n_r,plasma_sim.D_ambi(:,i),q);
    
    % Generate electron temperature diffusion PDE matrix, for fixed walls
    if plasma_sim.fixed_boundaries
        [M_t, bc_t] = matrix_T_input(n_r,plasma_sim.density(:,i),plasma_sim.k_therm(:,i),q,plasma_sim.Te_edge);
    end

    % Generate electron temperature diffusion PDE matrix, for sheath boundary (neumann bc)
    if plasma_sim.sheath_boundaries
        M_t = matrix_T_neumann(n_r,plasma_sim.density(:,i),plasma_sim.k_therm(:,i),q);
    end
    
    % Generate electron temperature diffusion PDE matrix, for heated boundary
    if plasma_sim.heated_boundaries
        ke_half = (1/2)*(plasma_sim.k_therm(end-1,i)+plasma_sim.k_therm(end,i));
        e_flux = ke_half*(plasma_sim.T_eles(end-1,i)-plasma_sim.T_eles(end,i))/dr;
        ki_half = (1/2)*(plasma_sim.ki_therm(end-1,i)+plasma_sim.ki_therm(end,i));
        i_flux = ki_half*(plasma_sim.T_ions(end-1,i)-plasma_sim.T_ions(end,i))/dr;
        plasma_sim.Te_edge = plasma_sim.Te_edge + (e_flux+i_flux)*dt/denom;
        plasma_sim.Ti_edge = plasma_sim.Te_edge;

        [M_t, bc_t] = matrix_T_input(n_r,plasma_sim.density(:,i),plasma_sim.k_therm(:,i),q,plasma_sim.Te_edge);
        
    end
    
    % Generate ion temperature diffusion PDE matrix
    [M_i,bc_i] = matrix_T_input(n_r,plasma_sim.density(:,i),plasma_sim.ki_therm(:,i),q,plasma_sim.Ti_edge);
    
end
toc;