% MATLAB script to compute Im[χ_0(q,ω)] in a 2D free electron gas
% Governs electr on-hole pair decay rates
clear; clc;

%% Define Physical Constants (in atomic units where ħ = m_e = 1)
Ef = 1.0;        % Fermi energy (sets energy scale) ~1 eV for typical metals
kF = sqrt(2*Ef); % Fermi wavevector (2D free electron gas)  
nq = 100;        % Number of q points
nw = 100;        % Number of frequency points

% Define q and ω ranges (normalized to kF and Ef)
q_vals = linspace(0, 2*kF, nq);   % Momentum transfer q (0 to 2*kF)
omega_vals = linspace(0, 2*Ef, nw); % Energy transfer ω (0 to 2*Ef)

% Initialize χ_0(q,ω)
ImChi0 = zeros(nq, nw);

%% Compute the Imaginary Part of χ_0(q,ω)
for iq = 1:nq
    q = q_vals(iq);
    for iw = 1:nw
        omega = omega_vals(iw);
        
        % Limits for electron-hole excitations (Landau damping region)
        q2 = q^2 / 2; % (q^2 / 2m, with m=1 in atomic units)
        omega_min = max(0, Ef - q2); % Lower bound for particle-hole excitations
        omega_max = Ef + q2;         % Upper bound
        
        % Im[χ_0(q,ω)] exists only in this window due to energy conservation
        if omega_min < omega && omega < omega_max
            ImChi0(iq, iw) = pi * Ef / (q * sqrt(omega_max^2 - omega^2));  
        else
            ImChi0(iq, iw) = 0; % No available states outside the range
        end
    end
end



%% Visualization of Im[χ_0(q,ω)] with Better Contrast
figure;
imagesc(q_vals/kF, omega_vals/Ef, log1p(ImChi0')); % log1p for better contrast
axis xy; 
colorbar;

% Use a high-contrast colormap (e.g., "plasma" or "magma")
colormap();

% Adjust color limits to enhance visibility
caxis([0 log1p(max(ImChi0(:))) * 0.8]); 

% Annotations
xlabel('q / k_F', 'FontSize', 12);
ylabel('\omega / E_F', 'FontSize', 12);
title('Imaginary Part of Electron Susceptibility Im[\chi_0(q,\omega)]', 'FontSize', 14);
set(gca, 'FontSize', 12);