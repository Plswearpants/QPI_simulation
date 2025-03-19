function [LDoS, worker_elapsed_time] = ComputeLDoS_vec( ...
    X_physical, omega, defect_energies, defect_locations_physical, ...
    a, t, E0, n, epsilon)
    % ComputeLDoS_vec
    %
    % This function computes a 2D Local Density of States (LDoS) on a grid
    % of observation points X_physical, given multiple energy values 'omega'
    % and a set of defects. The code uses vectorized operations to avoid
    % nested loops over defect pairs, and a parfor loop to parallelize
    % over the energy dimension. Each iteration is timed and the elapsed time 
    % is stored in 'worker_elapsed_time' and also printed.
    %
    % Inputs:
    %   X_physical: Grid of observation points (M x N x 2)
    %   omega:      Array of energy values of length num_energies
    %   defect_energies:          Array of defect energies
    %   defect_locations_physical: (num_defects x 2)
    %   a, t, E0, n, epsilon:     Physical parameters
    %
    % Outputs:
    %   LDoS:               (M x N x num_energies) local density of states
    %   worker_elapsed_time (num_energies x 1)     execution time per energy

    % Grid dimensions
    [M, N, ~]      = size(X_physical);
    num_defects    = size(defect_locations_physical, 1);
    num_energies   = length(omega);

    % Initialize output
    LDoS = zeros(M, N, num_energies);

    % Precompute T-matrix (defect-defect scattering) for all energies
    T_matrix = computeTMatrix(defect_energies, omega, ...
                              defect_locations_physical, a, t, E0, n, epsilon);

    % Optionally start a parallel pool outside or here
    if isempty(gcp('nocreate'))
        parpool('local');
    end

    % Store elapsed time for each energy
    worker_elapsed_time = zeros(num_energies,1);

    % Main parallel loop over energies
    parfor e = 1:num_energies
        % Start timing
        energy_tic = tic;

        % ------------------------------------------------------
        % 1) Compute G0_xd (grid->defects) and G0_dx (defects->grid)
        %    at this energy. 
        %    G0_xd: (M x N x num_defects)
        %    G0_dx: (num_defects x M x N), then we'll permute to (M x N x num_defects).
        % ------------------------------------------------------
        G0_xd = computeBLGF(X_physical, ...
                            defect_locations_physical, ...
                            omega(e), a, t, E0, n, epsilon);
        
        %G0_dx = computeBLGF(defect_locations_physical, ...
        %                    X_physical, ...
        %                    omega(e), a, t, E0, n, epsilon);
        %G0_dx = permute(G0_dx, [2, 3, 1]); % (M x N x num_defects)
        G0_dx = G0_xd;

        % ------------------------------------------------------
        % 2) Multiply G0_xd by T_matrix along the defect dimension
        %    T_ij is num_defects x num_defects for this energy e.
        % ------------------------------------------------------
        T_ij = T_matrix(:,:, e);   % (num_defects x num_defects)

        %   (M x N x num_defects) -> reshape -> (M*N) x num_defects
        g0_xd_reshaped = reshape(G0_xd, [M*N, num_defects]);

        %   multiply by T_ij -> (M*N) x num_defects
        g0_xd_times_T = g0_xd_reshaped * T_ij;

        %   reshape back -> (M x N x num_defects)
        g0_xd_times_T = reshape(g0_xd_times_T, [M, N, num_defects]);

        % ------------------------------------------------------
        % 3) Element-wise multiply by G0_dx and sum over defect dimension
        % ------------------------------------------------------
        worker_pre = sum( g0_xd_times_T .* G0_dx, 3 ); % (M x N)

        % ------------------------------------------------------
        % 4) Compute LDoS as (-1/pi)*Imag{worker_pre}
        % ------------------------------------------------------
        LDoS(:,:, e) = (-1/pi) * imag(worker_pre);

        % ------------------------------------------------------
        % Store and print timing
        % ------------------------------------------------------
        worker_elapsed_time(e) = toc(energy_tic);
        fprintf('Energy slice %d/%d (omega = %.4f): Computed in %.2f seconds\n', ...
            e, num_energies, omega(e), worker_elapsed_time(e));
    end
end
