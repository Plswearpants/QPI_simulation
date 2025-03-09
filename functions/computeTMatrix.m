function T_matrix = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon)
    % Compute T-matrix for defects in a lattice
    %
    % Inputs:
    %   defect_energies: Vector of defect energies [num_defects x 1]
    %   omega: Vector of energy values to compute
    %   defect_locations_physical: Defect positions [num_defects x 2]
    %   a: Lattice parameter
    %   t: Hopping parameter
    %   E0: Onsite energy
    %   n: Number of angle increments in Isq
    %   epsilon: Energy broadening
    % Outputs:
    %   T_matrix: 3D array [num_defects x num_defects x length(omega)]
    %             T_matrix(alpha,beta,i) gives the T-matrix element between
    %             defects alpha and beta at energy omega(i)
    num_defects = size(defect_locations_physical, 1);
    assert(length(defect_energies) == num_defects, 'Number of defect energies must match number of defects');
    
    % Compute G0 for each energy value
    % Note: defect_locations is now passed directly as Lx2 array
    G0 = arrayfun(@(i) computeBLGF(defect_locations_physical, defect_locations_physical, omega(i), a, t, E0, n, epsilon), ...
                  1:length(omega), 'UniformOutput', false);
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
        T_matrix(:,:,i) = E_alpha * (inner_term \ eye(num_defects));  % More numerically stable
    end
end
