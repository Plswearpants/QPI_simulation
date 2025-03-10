%% 1. Define Parameters
a = 1e-9;            % lattice constant (1 nm)
t = -0.2;           % hopping parameter
E0 = 0;             % on-site energy
n = 300;            % number of grid points for numerical integration
epsilon = 1e-3;     % small imaginary part for numerical stability
gridSize = 51;      % number of sampling points along one dimension
omega_values = linspace(-1, 1, 11);  % energy levels

%% 2. Single defect setup
defect_location = [0, 0];    % Single defect at origin
defect_energy = 0.002;         % Defect energy E1

%% 3. Create spatial grid for LDOS visualization
x_range = linspace(-10*a, 10*a, gridSize);
y_range = linspace(-10*a, 10*a, gridSize);
[X, Y] = meshgrid(x_range, y_range);
X_physical = cat(3, X, Y);
%% 4. Initialize arrays for results
T_single = zeros(1, length(omega_values));
delta_rho_1 = zeros(gridSize, gridSize, length(omega_values));
delta_rho_2 = zeros(gridSize, gridSize, length(omega_values));
delta_rho_3 = zeros(gridSize, gridSize, length(omega_values));
delta_rho_4 = zeros(gridSize, gridSize, length(omega_values));
delta_rho_5 = zeros(gridSize, gridSize, length(omega_values));

%% 5.1 G0_xd * T_single(i)* G0_dx
for i = 1:length(omega_values)
    % Compute G0(xd,xd) - Green's function at defect location
    G0_dd = computeBLGF(defect_location, defect_location, omega_values(i), a, t, E0, n, epsilon);
    
    % Compute T-matrix for single impurity
    T_single(i) = 1 / (1/defect_energy - G0_dd);
    
    % Compute G0(x,xd) for all grid points
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = squeeze(G0_dx);
    % Compute LDOS for all points
    delta_rho_1(:,:,i) = -1/pi * imag(G0_xd * T_single(i)* G0_dx);
end

%% 5.2 G0_xd * T_single(i)* G0_dx'
for i = 1:length(omega_values)
    % Compute G0(xd,xd) - Green's function at defect location
    G0_dd = computeBLGF(defect_location, defect_location, omega_values(i), a, t, E0, n, epsilon);
    
    % Compute T-matrix for single impurity
    T_single(i) = 1 / (1/defect_energy - G0_dd);
    
    % Compute G0(x,xd) for all grid points
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = squeeze(G0_dx);
    % Compute LDOS for all points
    delta_rho_2(:,:,i) = -1/pi * imag(G0_xd * T_single(i)* G0_dx');
end

% 5.3 G0_xd * T_single(i)* conj(G0_dx')
for i = 1:length(omega_values) 
    % Compute G0(xd,xd) - Green's function at defect location
    G0_dd = computeBLGF(defect_location, defect_location, omega_values(i), a, t, E0, n, epsilon);
    
    % Compute T-matrix for single impurity
    T_single(i) = 1 / (1/defect_energy - G0_dd);
    
    % Compute G0(x,xd) for all grid points
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = squeeze(G0_dx);
    % Compute LDOS for all points
    delta_rho_3(:,:,i) = -1/pi * imag(G0_xd * T_single(i)* conj(G0_dx'));
end

%% 5.4 G0_xd.*G0_xd * T_single(i)
for i = 1:length(omega_values)
    % Compute G0(xd,xd) - Green's function at defect location
    G0_dd = computeBLGF(defect_location, defect_location, omega_values(i), a, t, E0, n, epsilon);
    
    % Compute T-matrix for single impurity
    T_single(i) = 1 / (1/defect_energy - G0_dd);
    
    % Compute G0(x,xd) for all grid points
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = squeeze(G0_dx);
    % Compute LDOS for all points
    delta_rho_4(:,:,i) = -1/pi * imag(G0_xd.*G0_xd * T_single(i));
end

%% 5.5 G0_xd.*G0_xd' * T_single(i)
for i = 1:length(omega_values)
    % Compute G0(xd,xd) - Green's function at defect location
    G0_dd = computeBLGF(defect_location, defect_location, omega_values(i), a, t, E0, n, epsilon);
    
    % Compute T-matrix for single impurity
    T_single(i) = 1 / (1/defect_energy - G0_dd);
    
    % Compute G0(x,xd) for all grid points
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(i), a, t, E0, n, epsilon);
    G0_dx = squeeze(G0_dx);
    % Compute LDOS for all points
    delta_rho_5(:,:,i) = -1/pi * imag(G0_xd.*G0_xd' * T_single(i));
end

%% 6. Save all results
save('singledefect_results.mat', 'delta_rho_1', 'delta_rho_2', 'delta_rho_3', 'delta_rho_4', 'delta_rho_5');

