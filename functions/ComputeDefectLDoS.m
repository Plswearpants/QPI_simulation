function LDoS = ComputeDefectLDoS(X_physical, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
    [grid_x, grid_y] = size(X_physical(:,:,1));
    num_defects = size(defect_locations_physical, 1);
    
    % Compute G0 at each point
    G0 = arrayfun(@(i) computeBLGF(X_physical, defect_locations_physical, omega(i), a, t, E0, n, epsilon), 1:length(omega), 'UniformOutput', false);
    G0 = cat(4, G0{:});
    
    % Compute the multi-defect T-matrix
    T_matrix = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);
    
    % Compute the LDoS
    LDoS = zeros(grid_x, grid_y, length(omega));
    total_time = tic;
    for i = 1:length(omega)
        iter_time = tic;
        omega_val = omega(i);
        disp(['Computing LDoS for omega = ', num2str(omega_val)]);
        G0_i = G0(:,:,:,i);
        
        % Initialize sum for this frequency
        T_sum = zeros(grid_x, grid_y);
        
        % Sum over all defect pairs (alpha, beta)
        for alpha = 1:num_defects
            for beta = 1:num_defects
                % G0_i(:,:,alpha) gives G0(x,x_alpha) for all grid points x
                % G0_i(:,:,beta) gives G0(x_beta,x) for all grid points x
                T_element = T_matrix(alpha,beta,i);
                T_sum = T_sum + G0_i(:,:,alpha)* T_element* G0_i(:,:,beta);
            end
        end
        
        LDoS(:,:,i) = -imag(T_sum) / pi;
        elapsed = toc(iter_time);
        disp(['Iteration completed in ', num2str(elapsed), ' seconds']);
    end
    total_elapsed = toc(total_time);
    disp(['Total computation completed in ', num2str(total_elapsed), ' seconds']);
end
