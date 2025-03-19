function [convolved_LDoS] = compareDefectLDoS(single_defect_LDoS, defect_locations, N_multi, N_single, grid_size)
% Computes convolved LDoS from single defect pattern and multiple defect locations
%
% Inputs:
%   single_defect_LDoS - Grid of single defect LDoS (centered on 19x19 lattice)
%   defect_locations   - nx2 array of defect locations in lattice coordinates
%   N_lattice         - Size of the larger lattice (e.g., 50 for 50x50 lattice)
%   grid_size         - Size of the target sampling grid
%
% Output:
%   convolved_LDoS    - Resulting convolved LDoS on the grid_size x grid_size grid

    % Parameters for single defect pattern
    single_grid_size = size(single_defect_LDoS, 1);  % Get grid size from input array
    
    % Calculate scaling factors
    scale_factor_large = (grid_size - 1) / (N_multi - 1);  % For larger lattice
    scale_factor_single = (single_grid_size - 1) / (N_single - 1);  % For single defect
    
    % Calculate the physical size ratio between lattices
    physical_scale = scale_factor_single / scale_factor_large;
    
    % Resize single defect pattern to match physical size in larger lattice
    [X, Y] = meshgrid(1:single_grid_size, 1:single_grid_size);
    [Xq, Yq] = meshgrid(...
        linspace(1, single_grid_size, round(single_grid_size/physical_scale)), ...
        linspace(1, single_grid_size, round(single_grid_size/physical_scale)));
    single_defect_scaled = interp2(X, Y, single_defect_LDoS, Xq, Yq, 'cubic');
    
    % Get the size of the scaled pattern
    [single_height, single_width] = size(single_defect_scaled);
    
    % Initialize output array
    convolved_LDoS = zeros(grid_size, grid_size);
    
    % Convert lattice coordinates to grid coordinates
    grid_locations = (defect_locations - 1) * scale_factor_large + 1;
    
    % For each defect
    for i = 1:size(defect_locations, 1)
        % Get grid coordinates for this defect
        x_grid = round(grid_locations(i, 1));
        y_grid = round(grid_locations(i, 2));
        
        % Calculate the region where to place the scaled single defect pattern
        x_start = x_grid - floor(single_width/2);
        x_end = x_start + single_width - 1;
        y_start = y_grid - floor(single_height/2);
        y_end = y_start + single_height - 1;
        
        % Handle boundary conditions
        [x_start, x_end, y_start, y_end, pattern_x_start, pattern_x_end, ...
         pattern_y_start, pattern_y_end] = adjustBoundaries(x_start, x_end, ...
         y_start, y_end, grid_size, single_width, single_height);
        
        % Add the contribution of this defect
        convolved_LDoS(y_start:y_end, x_start:x_end) = ...
            convolved_LDoS(y_start:y_end, x_start:x_end) + ...
            single_defect_scaled(pattern_y_start:pattern_y_end, ...
                               pattern_x_start:pattern_x_end);
    end
    
    % Normalize the result
    %convolved_LDoS = convolved_LDoS / sum(convolved_LDoS(:));
end

function [x_start, x_end, y_start, y_end, pattern_x_start, pattern_x_end, ...
         pattern_y_start, pattern_y_end] = adjustBoundaries(x_start, x_end, ...
         y_start, y_end, grid_size, pattern_width, pattern_height)
    % Handle boundaries and calculate corresponding pattern indices
    
    pattern_x_start = 1;
    pattern_x_end = pattern_width;
    pattern_y_start = 1;
    pattern_y_end = pattern_height;
    
    % Adjust x boundaries
    if x_start < 1
        pattern_x_start = 2 - x_start;
        x_start = 1;
    end
    if x_end > grid_size
        pattern_x_end = pattern_width - (x_end - grid_size);
        x_end = grid_size;
    end
    
    % Adjust y boundaries
    if y_start < 1
        pattern_y_start = 2 - y_start;
        y_start = 1;
    end
    if y_end > grid_size
        pattern_y_end = pattern_height - (y_end - grid_size);
        y_end = grid_size;
    end
end 