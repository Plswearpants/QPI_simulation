function [LDoS_result,defect_locations] = computeLDoSWithMultipleDefects(a, t, E0, Ed, n, epsilon, num_defects, N, gridSize, omega_values)
    % Assign random defect locations without overlap
    defect_locations = assignDefectLocations(num_defects, N);

    % Scale defect locations to physical coordinates
    defect_locations_physical = defect_locations * a;

    % Initialize LDoS result
    LDoS_result = zeros(gridSize, gridSize, length(omega_values));

    % Create grid for physical coordinates
    [X, Y] = meshgrid(linspace(0, N*a, gridSize), linspace(0, N*a, gridSize));
    X = cat(3, X, Y);

    % Loop over omega values
    for i = 1:length(omega_values)
        omega = omega_values(i);
        tic;
        disp(['Computing LDoS for omega = ', num2str(omega)]);
        % Initialize G0 and T_matrix
        G0 = zeros(num_defects, num_defects);
        T_matrix = zeros(num_defects, num_defects);

        % Compute G0 and T_matrix for each pair of defect locations and
        % G0_xa_x and G0_x_xb
        for alpha = 1:num_defects
            for beta = alpha:num_defects
                if alpha == beta
                    G0(alpha, beta) = computeBLGF(zeros([1 1 2]), omega, a, t, E0, n, epsilon);
                    T_matrix(alpha, beta) = Ed / (1 - Ed * G0(alpha, beta));
                else
                    S = zeros([1 1 2]);
                    S(1,1,1) = (defect_locations_physical(alpha, 1) - defect_locations_physical(beta, 1));
                    S(1,1,2) = (defect_locations_physical(alpha, 2) - defect_locations_physical(beta, 2));
                    G0(alpha, beta) = computeBLGF(S, omega, a, t, E0, n, epsilon);
                    G0(beta, alpha) = G0(alpha, beta); % Symmetric property
                    T_matrix(alpha, beta) = -Ed * G0(alpha, beta) / (Ed * G0(beta, beta) - omega);
                    T_matrix(beta, alpha) = T_matrix(alpha, beta); % Symmetric property
                end
            end
        end

        % Compute G0 at each point
        [grid_x, grid_y] = size(X(:,:,1));
        S = zeros(grid_x, grid_y, 2, num_defects);
        for alpha = 1:num_defects
            S(:,:,1,alpha) = (X(:,:,1) - defect_locations_physical(alpha,1)) / a;
            S(:,:,2,alpha) = (X(:,:,2) - defect_locations_physical(alpha,2)) / a;
        end
        
        G0_points = zeros(grid_x, grid_y, num_defects);
        for alpha = 1:num_defects
            G0_points(:,:,alpha) = computeBLGF(S(:,:,:,alpha), omega, a, t, E0, n, epsilon);
        end
        % Compute the LDoS
        LDoS = zeros(grid_x, grid_y);
        for alpha = 1:num_defects
            for beta = alpha:num_defects
                if alpha == beta
                    LDoS = LDoS + G0_points(:,:,alpha) .* T_matrix(alpha, beta) .* G0_points(:,:,alpha);
                else
                    term = G0_points(:,:,alpha) .* T_matrix(alpha, beta) .* G0_points(:,:,beta);
                    LDoS = LDoS + 2* term ; % Add symmetric contribution
                end
            end
        end
        LDoS_result(:,:,i) = -imag(LDoS) / pi;
        toc;
    end
end



