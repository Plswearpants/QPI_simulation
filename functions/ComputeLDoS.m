function LDoS = ComputeLDoS(X_physical, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
    % This function computes a 2D Local density of state with defects
    %
    % Inputs:
    %   X_physical: Grid of observation points (MxNx2)
    %   omega: Array of energy values
    %   defect_energies: Array of defect energies
    %   defect_locations_physical: Array of defect locations (num_defects x 2)
    %   a, t, E0, n, epsilon: Physical parameters
    %
    % Output:
    %   LDoS: Local density of states (MxNxlength(omega))
    
    % Get grid dimensions
    [M, N, ~] = size(X_physical);
    
    % Initialize LDoS array
    LDoS = zeros(M, N, length(omega));
    
    % Compute T-matrix for all defects and energies
    T_matrix = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);
    
    % Compute LDoS for each energy value
    for e = 1:length(omega)
        tic; % Start timing for this energy slice
        
        % Compute Green's functions between grid points and defects
        G0_xd = computeBLGF(X_physical, defect_locations_physical, omega(e), a, t, E0, n, epsilon);
        G0_dx = computeBLGF(defect_locations_physical, X_physical, omega(e), a, t, E0, n, epsilon);
        
        % Initialize accumulator for this energy
        pre = zeros(M, N);
        
        % Sum contributions from all defect pairs
        for i = 1:size(defect_locations_physical, 1)
            for j = 1:size(defect_locations_physical, 1)
                % G0_xd(:,:,i) is MxN for the i-th defect
                % G0_dx(j,:,:) is 1xMxN for the j-th defect and needs to be squeezed
                pre = pre + G0_xd(:,:,i) .* squeeze(G0_dx(j,:,:)) * T_matrix(i,j,e);
            end
        end
        
        % Calculate LDOS change
        LDoS(:,:,e) = (-1/pi) * imag(pre);
        
        % Print the time taken for this energy slice
        elapsed_time = toc;
        fprintf('Energy slice %d/%d (ω = %.4f): Computed in %.2f seconds\n', e, length(omega), omega(e), elapsed_time);
    end
end
