function T_matrix = computeMultiDefectTMatrix(defect_energies, omega, defect_locations, a, t, E0, n, epsilon)
    num_defects = size(defect_locations, 1);
    
    % Compute G0 for each pair of defect locations
    [alpha, beta] = ndgrid(1:num_defects, 1:num_defects);
    dx = defect_locations(alpha, 1) - defect_locations(beta, 1);
    dy = defect_locations(alpha, 2) - defect_locations(beta, 2);
    S = cat(3, dx / a, dy / a);
    
    G0 = arrayfun(@(i) computeBLGF(S, omega(i), a, t, E0, n, epsilon), 1:length(omega), 'UniformOutput', false);
    G0 = cat(3, G0{:});
    
    % Compute the T-matrix
    T_matrix = zeros(num_defects, num_defects, length(omega));
    for i = 1:length(omega)
        omega_i = omega(i);
        G0_i = G0(:,:,i);
        T_matrix(:,:,i) = diag(defect_energies ./ (1 - defect_energies .* diag(G0_i))) - ...
                          defect_energies .* G0_i ./ (defect_energies .* diag(G0_i) - omega_i);
    end
end
