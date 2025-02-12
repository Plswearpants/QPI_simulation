function [LDoS_result,defect_locations] = test_multildos(a, t, E0, Ed, n, epsilon, num_defects, N, gridSize, omega_values)
    % Assign random defect locations without overlap
    defect_locations = assignDefectLocations(num_defects, N);
    
    % Scale defect locations to physical coordinates
    defect_locations_physical = defect_locations * a;
    
    % Initialize LDoS result
    LDoS_result = zeros(gridSize, gridSize, length(omega_values));
    
    % Create grid for physical coordinates
    [X, Y] = meshgrid(linspace(0, N*a, gridSize), linspace(0, N*a, gridSize));
    
    % Loop over omega values
    for i = 1:length(omega_values)
        omega = omega_values(i);
        tic;
        disp(['Computing LDoS for omega = ', num2str(omega)]);
        
        % Step 1: Compute G₀(xα,xβ) and T-matrix
        G0 = zeros(num_defects, num_defects);
        T_matrix = zeros(num_defects, num_defects);
        
        % Compute G₀(xα,xβ) for all pairs
        for alpha = 1:num_defects
            for beta = 1:num_defects
                if alpha == beta
                    G0(alpha, beta) = computeBLGF(zeros([1 1 2]), omega, a, t, E0, n, epsilon);
                else
                    S = zeros([1 1 2]);
                    S(1,1,1) = (defect_locations_physical(alpha, 1) - defect_locations_physical(beta, 1));
                    S(1,1,2) = (defect_locations_physical(alpha, 2) - defect_locations_physical(beta, 2));
                    G0(alpha, beta) = computeBLGF(S, omega, a, t, E0, n, epsilon);
                end
            end
        end
        
        % Compute T-matrix
        % T = V(1-G₀V)^(-1), where V is the defect potential
        V = Ed * eye(num_defects);
        T_matrix = V * inv(eye(num_defects) - G0 * V);
        
        % Step 2: Compute G₀(x,xα) for all grid points and defects
        G0_x_xa = zeros(gridSize, gridSize, num_defects);
        for alpha = 1:num_defects
            S = zeros(gridSize, gridSize, 2);
            S(:,:,1) = X - defect_locations_physical(alpha,1);
            S(:,:,2) = Y - defect_locations_physical(alpha,2);
            G0_x_xa(:,:,alpha) = computeBLGF(S/a, omega, a, t, E0, n, epsilon);
        end
        
        % Step 3: Compute LDoS according to the equation
        LDoS = zeros(gridSize, gridSize);
        for alpha = 1:num_defects
            for beta = 1:num_defects
                LDoS = LDoS + G0_x_xa(:,:,alpha) .* T_matrix(alpha,beta) .* conj(G0_x_xa(:,:,beta));
            end
        end
        
        % Final step: Apply the formula δρ(x,ω) = -(1/π) Im[...]
        LDoS_result(:,:,i) = -imag(LDoS) / pi;
        toc;
    end
end