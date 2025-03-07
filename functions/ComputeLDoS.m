function LDoS = ComputeLDoS(X_physical, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
    % This function computes a 2D Local density of state with single defect.
    
    % Empty LDoS
    LDoS = zeros(size(X_physical, 1), size(X_physical, 2), length(omega));
    
    % Compute T-matrix
    T = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);
    % Compute LDoS 
    for i = 1:length(omega)
        G0 = computeBLGF(X_physical, defect_locations_physical, omega(i), a, t, E0, n, epsilon);
        LDoS(:,:,i) = -imag(G0*G0* T(1,1,i)) / pi;
    end
end
