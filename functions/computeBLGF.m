function G0 = computeBLGF(S, omega, a, t, E0, n, epsilon)
    
arguments
    S       % seperation vector between 2D vectors X and X', S= X-X'
    omega   % energy slice to compute
    a       % lattice parameter
    t       % hopping parameters    
    E0      % onsite energy of the lattice 
    n       % number of angle incrementation in Isq. larger n, more accurate Isq integration 
    epsilon % broadening of the energy
end

    % Bare Lattice Green's Function Calculation, outputs a 2D G0 values
    s1 = S(:,:,1);
    s2 = S(:,:,2);
    G0 = zeros(size(S, 1), size(S, 2));
    G0_flat = zeros(numel(s1), 1);
    
    % Vectorize the integration over the grid
    b = (omega + 1i*epsilon+1i*0.0627*omega^2/2 - E0) / (2 * t) ; % Add small imaginary part for numerical stability
    
    % Flatten s1 and s2 for parallel computing 
    s1_flat = s1(:);
    s2_flat = s2(:);
    
    % Compute the integral for each point using arrayfun
    G0_flat(:) = arrayfun(@(bx, by) Isq_integral(bx, by, b, n), s1_flat, s2_flat);
    
    % Normalize and reshape back to 2D grid
    G0_flat = G0_flat * 4 / ((2 * pi * a)^2 * 2 * t);
    G0 = reshape(G0_flat, size(S, 1), size(S, 2));
end