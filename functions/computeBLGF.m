function G0 = computeBLGF(X1, X2, omega, a, t, E0, n, epsilon)
    
arguments
    X1      % First position vector (x,y)
    X2      % Second position vector (x',y')
    omega   % energy slice to compute
    a       % lattice parameter
    t       % hopping parameters    
    E0      % onsite energy of the lattice 
    n       % number of angle incrementation in Isq. larger n, more accurate Isq integration 
    epsilon % broadening of the energy
end

    % Compute separation vector S = (X1-X2)/a according to the formula
    S = zeros(size(X1));
    S(:,:,1) = (X1(:,:,1) - X2(:,:,1)) / a;  % s₁ = (x₁ - x'₁)/a
    S(:,:,2) = (X1(:,:,2) - X2(:,:,2)) / a;  % s₂ = (x₂ - x'₂)/a
    
    % Bare Lattice Green's Function Calculation
    s1 = S(:,:,1);
    s2 = S(:,:,2);
    G0 = zeros(size(S, 1), size(S, 2));
    G0_flat = zeros(numel(s1), 1);
    
    % Vectorize the integration over the grid
    b = (omega + 1i*epsilon+1i*0.0627*omega^2/2 - E0) / (2 * t);
    
    % Flatten s1 and s2 for parallel computing 
    s1_flat = s1(:);
    s2_flat = s2(:);
    
    % Compute the integral for each point using arrayfun
    G0_flat(:) = arrayfun(@(bx, by) Isq_integral(bx, by, b, n), s1_flat, s2_flat);
    
    % Normalize and reshape back to 2D grid
    G0_flat = G0_flat * 4 / ((2 * pi * a)^2 * 2 * t);
    G0 = reshape(G0_flat, size(S, 1), size(S, 2));
end