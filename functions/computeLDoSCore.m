function [LDoS, defect_locations_lattice] = computeLDoSCore(omega, defect_energies, defect_locations_lattice, N, a, t, E0, n, epsilon, gridSize)
    % computeLDoSCore - Computes LDoS for given defect locations
    %
    % Part 1: Input validation and physical grid creation
    % Part 2: Computation of LDoS using ComputeDefectLDoS
    %
    % Inputs:
    %   omega - Array of frequency values
    %   defect_energies - Array of defect energies
    %   defect_locations - Matrix of defect locations (Nx2): location in the lattice, not physical position.
    %   N - Number of lattice points along one dimension
    %   a - Lattice constant
    %   t - Hopping parameter
    %   E0 - Energy parameter
    %   n - Integration parameter
    %   epsilon - Small parameter for integration
    %   gridSize - Size of the computation grid
    %
    % Outputs:
    %   LDoS - Local Density of States
    %   defect_locations - Used defect locations
    
    % Part 1: Validate inputs and create physical grid
    if size(defect_locations_lattice, 2) ~= 2
        error('defect_locations must be a Nx2 matrix');
    end
    if length(defect_energies) ~= size(defect_locations_lattice, 1)
        error('Number of defect energies must match number of defect locations');
    end

    % Create physical coordinate grid
    [X1, X2] = meshgrid(linspace(-N*a/2, N*a/2, gridSize), linspace(-N*a/2, N*a/2, gridSize));
    X = cat(3, X1, X2); % Location vector on the grid

    % Convert defect locations from lattice coordinates to physical coordinates
    % Shift from (N/2 + 0.5, N/2 + 0.5) centered coordinates to (0,0) centered coordinates
    defect_locations_physical = (defect_locations_lattice - (N/2 + 0.5)) * a;

    % Part 2: Compute LDoS using the core computation function
    %LDoS = ComputeDefectLDoS(X, omega, defect_energies, defect_locations, ...
    %                       a, t, E0, n, epsilon);
    LDoS = ComputeLDoS(X, omega, defect_energies, defect_locations_physical, a, t, E0, n, epsilon);
end 