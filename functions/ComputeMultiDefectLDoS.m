function LDoS = ComputeMultiDefectLDoS(X, omega, defect_energies, defect_locations, a, t, E0, n, epsilon)
    [grid_x, grid_y] = size(X(:,:,1));
    num_defects = size(defect_locations, 1);
    
    % Compute G0 at each point
    [x, y] = ndgrid(1:grid_x, 1:grid_y);
    S = zeros(grid_x, grid_y, 2, num_defects);
    S(:,:,1,:) = reshape((x(:) - defect_locations(:,1)') / a, [grid_x, grid_y, 1, num_defects]);
    S(:,:,2,:) = reshape((y(:) - defect_locations(:,2)') / a, [grid_x, grid_y, 1, num_defects]);
    
    G0 = arrayfun(@(i) computeBLGF(S, omega(i), a, t, E0, n, epsilon), 1:length(omega), 'UniformOutput', false);
    G0 = cat(4, G0{:});
    
    % Compute the multi-defect T-matrix
    T_matrix = computeMultiDefectTMatrix(defect_energies, omega, defect_locations, a, t, E0, n, epsilon);
    
    % Compute the LDoS
    tic;
    LDoS = zeros(grid_x, grid_y, length(omega));
    for i = 1:length(omega)
        omega_val = omega(i);
        disp(['Computing LDoS for omega = ', num2str(omega_val)]);
        G0_i = G0(:,:,:,i);
        T_sum = sum(sum(G0_i .* permute(T_matrix(:,:,i), [3, 4, 1, 2]) .* permute(G0_i, [1, 2, 4, 3]), 4), 3);
        LDoS(:,:,i) = -imag(T_sum) / pi;
        toc
    end
    elapsed_time = toc;
end
