function G0 = computeBLGF(R1, R2, omega, a, t, E0, n, epsilon)
    
arguments
    R1      % Position vectors in physical coordinates - either grid (MxNx2) or list of points (Kx2)
    R2      % Position vectors in physical coordinates - list of points (Lx2)
    omega   % Energy slice to compute
    a       % Lattice parameter
    t       % Hopping parameters    
    E0      % Onsite energy of the lattice 
    n       % Number of angle incrementation in Isq
    epsilon % Broadening of the energy
end
    % Output G0: Green's function between points in R1 and R2
    % In one word, dim(G0) = dim(R1) + dim(R2). 
    % Case 1: When R1 is a grid (MxNx2)
    %   G0 is a 3D array of size [M, N, L] where:
    %   - M: number of rows in R1 grid
    %   - N: number of columns in R1 grid
    %   - L: number of points in R2
    %   G0(i,j,k) gives the Green's function between grid point R1(i,j,:) 
    %   and defect point R2(k,:)
    %
    % Case 2: When R1 is a list (Kx2)
    %   G0 is a 2D array of size [K, L] where:
    %   - K: number of points in R1
    %   - L: number of points in R2
    %   G0(i,j) gives the Green's function between point R1(i,:) 
    %   and point R2(j,:)
    % Store original dimensions of R1 if it's a grid

    if ndims(R1) == 3  
        [M, N, ~] = size(R1);
        R1_list = reshape(R1, M*N, 2);  % Convert grid to list of points
    else
        R1_list = R1;
    end
    
    % Get number of points
    n_points_R1 = size(R1_list, 1);
    n_points_R2 = size(R2, 1);
    
    % Initialize separation vectors
    S = zeros(n_points_R1, n_points_R2, 2);
    
    % Compute separation vectors between all pairs of points
    S(:,:,1) = (R1_list(:,1) - R2(:,1)') / a;  % x components
    S(:,:,2) = (R1_list(:,2) - R2(:,2)') / a;  % y components
    
    % Extract components for Green's function calculation
    sx = S(:,:,1);
    sy = S(:,:,2);
    
    % Initialize Green's function arrays
    G0 = zeros(n_points_R1, n_points_R2);
    G0_flat = zeros(numel(sx), 1);
    
    % Compute Green's function
    b = (omega + 1i*epsilon + 1i*0.0627*omega^2/2 - E0) / (2 * t);
    
    % Flatten arrays for parallel computing
    sx_flat = sx(:);
    sy_flat = sy(:);
    
    % Compute integral for each point
    G0_flat(:) = arrayfun(@(sx, sy) Isq_integral(sx, sy, b, n), sx_flat, sy_flat);
    
    % Normalize and reshape
    G0_flat = G0_flat * 4 / ((2 * pi * a)^2 * 2 * t);
    G0 = reshape(G0_flat, n_points_R1, n_points_R2);
    
    % If R1 was originally a grid, reshape G0 back to grid form
    if ndims(R1) == 3
        G0 = reshape(G0, [M, N, n_points_R2]);  % Explicitly specify dimensions [M,N,L]
    end
end