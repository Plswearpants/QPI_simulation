function LDoS = ComputeLDoS(X_physical, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon)
    % This function computes a 2D Local density of state with defects
    
    % Get grid dimensions
    [M, N, ~] = size(X_physical);
    
    % Empty LDoS
    LDoS = zeros(M, N, length(omega));
    
    % Compute T-matrix
    T = computeTMatrix(defect_energies, omega, defect_locations_physical, a, t, E0, n, epsilon);
    
    % Compute LDoS 
    for i = 1:length(omega)
        % Compute G0(x,xₐ) - from observation points to defects
        G0_xd = computeBLGF(X_physical, defect_locations_physical, omega(i), a, t, E0, n, epsilon);

        % Compute G0(xᵦ,x) - from defects to observation points
        % Note: this is the conjugate transpose of G0(x,xᵦ) due to reciprocity
        G0_dx = G0_xd';
        
        % Compute full sum: G0(x,xₐ)·Tₐᵦ·G0(xᵦ,x)
        LDoS(:,:,i) = -imag(G0_xd * T(:,:,i) * G0_dx) / pi;
        size(LDoS(:,:,i))
    end
end
