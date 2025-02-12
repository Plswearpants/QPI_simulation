function LDoS = ComputeLDoS(X, omega, a, t, E0, Ed, Xd, n, epsilon)
    % This function computes a 2D Local density of state with single defect.
    
    % Displacement grid
    S = zeros(size(X));
    S(:,:,1) = (X(:,:,1) - Xd(1)) / a;
    S(:,:,2) = (X(:,:,2) - Xd(2)) / a;
    
    % Empty LDoS
    LDoS = zeros(size(X, 1), size(X, 2), length(omega));
    
    % Compute T-matrix
    T = computeTMatrix(Ed, omega, a, t, E0, n, epsilon);
   
    % Compute LDoS 
    for i = 1:length(omega)
        G0 = computeBLGF(S, omega(i), a, t, E0, n, epsilon);
        LDoS(:,:,i) = -imag(G0 .* G0 * T(i)) / pi;
    end
end
