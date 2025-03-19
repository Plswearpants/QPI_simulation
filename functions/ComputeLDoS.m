function [LDoS, worker_LDoS] = ComputeLDoS(X_physical, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
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
    num_defects = size(defect_locations_physical, 1);
    num_energies = length(omega);
    
    % Initialize LDoS array - this will collect results from all workers
    LDoS = zeros(M, N, num_energies);
    
    % Compute T-matrix for all defects and energies
    % This is done outside the parfor loop to avoid redundant computation
    T_matrix = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);

    % Initialize parallel pool if not already running
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    
    %pre-allocate memory for worker_LDoS, worker_pre, worker_G0_xd, worker_G0_dx, and worker_T_ij, worker_elapsed_time
    worker_LDoS = cell(num_energies, 1);
    worker_pre = cell(num_energies, 1);
    worker_G0_xd = cell(num_energies, 1);
    worker_G0_dx = cell(num_energies, 1);
    worker_T_ij = cell(num_energies, 1);
    worker_elapsed_time = cell(num_energies, 1);
    worker_timer = cell(num_energies, 1);
    % Use parallel computing for energy values
    parfor e = 1:num_energies
        % Each worker gets its own copy of these variables
        worker_pre{e} = zeros(M, N);
        worker_G0_xd{e} = zeros(M, N, num_defects);
        worker_G0_dx{e} = zeros(num_defects, M, N);
        worker_T_ij{e} = T_matrix(:,:,e);
        
        % Start timing for this energy slice
        worker_timer{e} = tic;
        
        % Compute Greens functions between grid points and defects
        % Each worker gets its own copy of these matrices
        worker_G0_xd{e} = computeBLGF(X_physical, defect_locations_physical, omega(e), a, t, E0, n, epsilon);
        worker_G0_dx{e} = computeBLGF(defect_locations_physical, X_physical, omega(e), a, t, E0, n, epsilon);
        
        % Initialize accumulator for this energy - each worker has its own
        worker_pre{e} = zeros(M, N);
        
        % Sum contributions from all defect pairs
        for i = 1:num_defects
            for j = 1:num_defects
                % Accumulate contribution to workers local accumulator
                worker_pre{e} = worker_pre{e} + worker_G0_xd{e}(:,:,i) .* squeeze(worker_G0_dx{e}(j,:,:)) * worker_T_ij{e}(i,j);
            end
        end
        
        % Calculate LDOS change - each worker computes its own slice
        LDoS(:,:,e) = (-1/pi) * imag(worker_pre{e});
        
        % Print the time taken for this energy slice
        worker_elapsed_time{e} = toc(worker_timer{e});
        fprintf('Energy slice %d/%d (ω = %.4f): Computed in %.2f seconds\n', ...
                e, num_energies, omega(e), worker_elapsed_time{e});
    end
end
