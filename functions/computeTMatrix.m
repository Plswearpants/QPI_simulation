function T = computeTMatrix(E, omega, a, t, E0, n, epsilon)
    % T-Matrix Calculation for Single Defect
    % Compute G0 at the defect position (0 displacement)
    G0_xdxd = computeBLGF(zeros(size(omega, 1), size(omega, 2), 2), omega, a, t, E0, n, epsilon);
    T = 1 ./ (E - G0_xdxd);
end