function LDoS = ComputeLDoS(X, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
    % This function computes a 2D Local density of state with single defect.
    
    % Displacement grid
    S = zeros(size(X));
    S(:,:,1) = (X(:,:,1) - defect_locations_physical(1)) / a;
    S(:,:,2) = (X(:,:,2) - defect_locations_physical(2)) / a;
    
    % Empty LDoS
    LDoS = zeros(size(X, 1), size(X, 2), length(omega));
    
    % Compute T-matrix
    T = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);
    % Compute LDoS 
    for i = 1:length(omega)
        G0 = computeBLGF(S, omega(i), a, t, E0, n, epsilon);
        LDoS(:,:,i) = -imag(G0 .* G0 * T(1,1,i)) / pi;
    end
end
