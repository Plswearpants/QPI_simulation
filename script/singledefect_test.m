%% 1. Define Parameters
a = 1e-9;            % lattice constant (1 nm)
t = -0.2;           % hopping parameter
E0 = 0;             % on-site energy
n = 300;            % number of grid points for numerical integration
epsilon = 1e-3;     % small imaginary part for numerical stability
gridSize = 128;      % number of sampling points along one dimension
omega_values = linspace(-0.5, 0.5, 21);  % energy levels

%% 2. Single defect setup
defect_location = [[0, 0];[1*a,2*a]];    % Single defect at origin
defect_energy = 0.002;         % Defect energy E1

%% 3. Create spatial grid for LDOS visualization
x_range = linspace(-15*a, 15*a, gridSize);
y_range = linspace(-15*a, 15*a, gridSize);
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
    delta_rho_4(:,:,i) = -1/pi * imag(G0_xd.*G0_dx * T_single(i));
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


%% 6. Multi-defect LDOS calculation
% Initialize array for multi-defect LDOS results
delta_rho_multi = zeros(gridSize, gridSize, length(omega_values));

% Define defect parameters
defect_energies = 0.002; % Energy for each defect
%defect_location = [[0, 0]; [4*a, 5*a]; [-10*a, -2*a]]; % Multiple defect locations
defect_location = [0, 0]; % single defect locations

% Start timing the entire calculation
total_time_start = tic;

% Compute T-matrix for multiple defects
T_matrix = computeTMatrix(defect_energies, omega_values, defect_location, a, t, E0, n, epsilon);

% Calculate LDOS for each energy value using original nested loop method
for e = 1:length(omega_values)
    % Start timing this energy slice
    energy_time_start = tic;
    
    % Compute Green's functions between grid points and defects
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(e), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(e), a, t, E0, n, epsilon);
    
    % Initialize accumulator for this energy
    pre = zeros(gridSize, gridSize);
    
    % Sum contributions from all defect pairs
    for i = 1:size(defect_location, 1)
        for j = 1:size(defect_location, 1)
            % Note: squeeze is important to handle dimensions correctly
            % G0_xd(:,:,i) is gridSize x gridSize
            % G0_dx(j,:,:) is 1 x gridSize x gridSize and needs to be squeezed
            pre = pre + G0_xd(:,:,i) .* squeeze(G0_dx(j,:,:)) * T_matrix(i,j,e);
        end
    end
    
    % Calculate LDOS change
    delta_rho_multi(:,:,e) = (-1/pi) * imag(pre);
    
    % End timing for this energy slice and display
    energy_time_elapsed = toc(energy_time_start);
    fprintf('Energy slice %d/%d completed in %.4f seconds\n', e, length(omega_values), energy_time_elapsed);
end

% End timing for the entire calculation and display
total_time_elapsed = toc(total_time_start);
fprintf('Total multi-defect LDOS calculation completed in %.4f seconds\n', total_time_elapsed);

%% 6.1 Optimized Multi-defect LDOS calculation
omega_values = linspace(-0.5,0.5,3);
defect_location = [[0, 0]; [4*a, 5*a]; [-10*a, -2*a]]; % Multiple defect locations
defect_energies = [0.002,0.002,0.002]; % Energy for each defect

% Initialize array for optimized multi-defect LDOS results
delta_rho_multi_opt = zeros(gridSize, gridSize, length(omega_values));

% Start timing the entire optimized calculation
opt_total_time_start = tic;

% Compute T-matrix for multiple defects
T_matrix = computeTMatrix(defect_energies, omega_values, defect_location, a, t, E0, n, epsilon);

% Calculate LDOS for each energy value using optimized vectorized method
for e = 1:length(omega_values)
    % Start timing this energy slice
    opt_energy_time_start = tic;
    
    % Compute Green's functions between grid points and defects
    G0_xd = computeBLGF(X_physical, defect_location, omega_values(e), a, t, E0, n, epsilon);
    G0_dx = computeBLGF(defect_location, X_physical, omega_values(e), a, t, E0, n, epsilon);
    
    % Vectorized approach using matrix operations
    num_defects = size(defect_location, 1);
    
    % Reshape G0_xd to prepare for matrix multiplication
    G0_xd_reshaped = reshape(G0_xd, [gridSize*gridSize, num_defects]);
    
    % Reshape G0_dx to prepare for element-wise multiplication
    G0_dx_reshaped = permute(G0_dx, [2, 3, 1]); % Rearrange to [gridSize, gridSize, num_defects]
    G0_dx_reshaped = reshape(G0_dx_reshaped, [gridSize*gridSize, num_defects]);
    
    % Compute the sum using matrix multiplication
    T_e = T_matrix(:,:,e); % Extract T-matrix for current energy
    result = G0_xd_reshaped * T_e .* G0_dx_reshaped;
    
    % Sum over all defect contributions and reshape back to grid
    pre_vectorized = reshape(sum(result, 2), [gridSize, gridSize]);
    
    % Calculate LDOS change
    delta_rho_multi_opt(:,:,e) = (-1/pi) * imag(pre_vectorized);
    
    % End timing for this energy slice and display
    opt_energy_time_elapsed = toc(opt_energy_time_start);
    fprintf('Optimized energy slice %d/%d completed in %.4f seconds\n', e, length(omega_values), opt_energy_time_elapsed);
end

% End timing for the entire optimized calculation and display
opt_total_time_elapsed = toc(opt_total_time_start);
fprintf('Total optimized multi-defect LDOS calculation completed in %.4f seconds\n', opt_total_time_elapsed);

%% 6.2 Verify optimization correctness
% Test 1: Compare the results from both methods
max_diff = max(abs(delta_rho_multi(:) - delta_rho_multi_opt(:)));
fprintf('Maximum difference between original and optimized methods: %e\n', max_diff);

% Test 2: Relative difference as percentage
rel_diff = 100 * max_diff / max(abs(delta_rho_multi(:)));
fprintf('Maximum relative difference: %.10f%%\n', rel_diff);

% Test 3: Check if results are effectively identical (within numerical precision)
is_identical = max_diff < 1e-10;
fprintf('Results are %s\n', conditional(is_identical, 'identical within numerical precision', 'different'));

% Helper function for conditional string output
function result = conditional(condition, true_str, false_str)
    if condition
        result = true_str;
    else
        result = false_str;
    end
end

%% 7. Save all results
%save('singledefect_results_1.mat', 'delta_rho_1', 'delta_rho_2', 'delta_rho_3', 'delta_rho_4', 'delta_rho_5', 'delta_rho_multi', 'delta_rho_multi_opt');
save('multidefect_results.mat','delta_rho_multi');

%% 8. Visualize multi-defect LDOS comparison (optional)
figure;
for i = 1:length(omega_values)
    % Original method
    subplot(2, length(omega_values), i);
    imagesc(x_range/a, y_range/a, delta_rho_multi(:,:,i));
    colorbar;
    title(['Original: \omega = ', num2str(omega_values(i))]);
    xlabel('x/a');
    ylabel('y/a');
    axis equal tight;
    
    % Optimized method
    subplot(2, length(omega_values), i + length(omega_values));
    imagesc(x_range/a, y_range/a, delta_rho_multi_opt(:,:,i));
    colorbar;
    title(['Optimized: \omega = ', num2str(omega_values(i))]);
    xlabel('x/a');
    ylabel('y/a');
    axis equal tight;
end

% Plot difference (if any)
if ~is_identical
    figure;
    for i = 1:length(omega_values)
        subplot(1, length(omega_values), i);
        imagesc(x_range/a, y_range/a, abs(delta_rho_multi(:,:,i) - delta_rho_multi_opt(:,:,i)));
        colorbar;
        title(['Difference: \omega = ', num2str(omega_values(i))]);
        xlabel('x/a');
        ylabel('y/a');
        axis equal tight;
    end
end

