% MATLAB script to compute Im[χ_RPA(q,ω)] and simulate temporal decay
clear; clc;

%% Define Physical Constants (atomic units where ħ = m_e = 1)
Ef = 1.0;         % Fermi energy (eV)
kF = sqrt(2*Ef);  % Fermi wavevector
e2 = 1.0;         % Effective Coulomb coupling (scaled for simplicity)
epsilon = 2.0;    % Background dielectric constant (screening)

nq = 100;  % Number of q points
nw = 100;  % Number of frequency points
nt = 100;  % Number of time points

% Define q and ω ranges
q_vals = linspace(0.1, 2*kF, nq);    % Avoid q=0 to prevent division issues
omega_vals = linspace(0, 2*Ef, nw);  % Energy range

% Initialize susceptibility matrices
ImChi0 = zeros(nq, nw);
ImChiRPA = zeros(nq, nw);
Gamma = zeros(nq, nw);

%% Compute Non-Interacting Im[χ_0(q,ω)]
for iq = 1:nq
    q = q_vals(iq);
    for iw = 1:nw
        omega = omega_vals(iw);
        
        % Limits for electron-hole excitations (Landau damping)
        q2 = q^2 / 2; 
        omega_min = max(0, Ef - q2); 
        omega_max = Ef + q2; 
        
        if omega_min < omega && omega < omega_max
            ImChi0(iq, iw) = pi * Ef / (q * sqrt(omega_max^2 - omega^2));  
        else
            ImChi0(iq, iw) = 0; 
        end
    end
end

%% Compute RPA-Corrected Im[χ_RPA(q,ω)]
for iq = 1:nq
    q = q_vals(iq);
    Vq = (2 * pi * e2) / (epsilon * q); % 2D Coulomb interaction
    for iw = 1:nw
        omega = omega_vals(iw);
        
        % RPA correction: χ_RPA = χ_0 / (1 - V(q) χ_0)
        denominator = 1 - Vq * ImChi0(iq, iw);
        if abs(denominator) > 1e-6
            ImChiRPA(iq, iw) = ImChi0(iq, iw) / denominator;
        else
            ImChiRPA(iq, iw) = 0; % Avoid singularities
        end
    end
end

% Compute decay rate Gamma(q, ω)
Gamma = -2 * ImChiRPA;

%% Simulate Temporal Decay S(q,t)
t_vals = linspace(0, 10, nt); % Time range
S_q_t = zeros(nq, nt);

for iq = 1:nq
    for it = 1:nt
        t = t_vals(it);
        % Fourier transform to get time evolution
        S_q_t(iq, it) = trapz(omega_vals, Gamma(iq, :) .* exp(-1i * omega_vals * t));
    end
end

%% Visualization with Fixed Complex Values
figure;
subplot(1,2,1);
imagesc(q_vals/kF, omega_vals/Ef, log1p(abs(ImChiRPA)')); % Ensuring real values
axis xy; colorbar; colormap;
xlabel('q / k_F'); ylabel('\omega / E_F');
title('RPA-Corrected Im[\chi(q,\omega)]');

subplot(1,2,2);
imagesc(q_vals/kF, t_vals, real(S_q_t)'); % Extract real part
axis xy; colorbar; colormap turbo;
xlabel('q / k_F'); ylabel('Time (arb. units)');
title('Temporal Decay of Electron-Hole Pair S(q,t)');
