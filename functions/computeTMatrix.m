function T_matrix = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon)
    num_defects = size(defect_locations_physical, 1);
    assert(length(defect_energies) == num_defects, 'Number of defect energies must match number of defects');
    
    % Compute G0 for each pair of defect locations
    [alpha, beta] = ndgrid(1:num_defects, 1:num_defects);
    
    G0 = arrayfun(@(i) computeBLGF(defect_locations_physical(alpha, :), defect_locations_physical(beta, :), omega(i), a, t, E0, n, epsilon), 1:length(omega), 'UniformOutput', false);
    G0 = cat(3, G0{:});
    
    % Compute the T-matrix
    T_matrix = zeros(num_defects, num_defects, length(omega));
    for i = 1:length(omega)
        G0_i = G0(:,:,i);
        
        % Create matrices for the calculation
        E_alpha = diag(defect_energies);  % E_α matrix
        E_beta = diag(defect_energies);   % E_β matrix
        
        % Calculate (δ_αβ - E_β G₀(x_α, x_β))
        delta = eye(num_defects);  % Kronecker delta δ_αβ
        inner_term = delta - E_beta * G0_i;
        
        % Calculate T_αβ = E_α (δ_αβ - E_β G₀(x_α, x_β))⁻¹
        T_matrix(:,:,i) = E_alpha * (inner_term \ eye(num_defects));
    end
end
