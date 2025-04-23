function [multi_defect_LDoS, defect_locations] = simulate_multi_defect(single_defect_LDoS, N, num_defects, varargin)
    % SIMULATE_MULTI_DEFECT Simulates multi-defect LDoS by convolving single defect LDoS
    % with random defect positions on a lattice
    %
    % Inputs:
    %   single_defect_LDoS - 2D array of single defect LDoS
    %   N - number of lattice sites in each dimension
    %   num_defects - number of defects to place
    %   varargin - Optional name-value pairs:
    %       'PlotResult' - logical flag to show visualization (default: true)
    %       'DefectLocations' - [x,y] coordinates for defects (if not random)
    %
    % Outputs:
    %   multi_defect_LDoS - simulated multi-defect LDoS
    %   defect_locations - [x,y] coordinates of placed defects
    
    % Parse inputs
    p = inputParser;
    addRequired(p, 'single_defect_LDoS', @isnumeric);
    addRequired(p, 'N', @isnumeric);
    addRequired(p, 'num_defects', @isnumeric);
    addParameter(p, 'PlotResult', true, @islogical);
    addParameter(p, 'DefectLocations', [], @isnumeric);
    parse(p, single_defect_LDoS, N, num_defects, varargin{:});
    
    % Get grid size
    grid_size = size(single_defect_LDoS, 1);
    scale_factor = (grid_size - 1) / (N - 1);
    
    % Generate random defect locations if not provided
    if isempty(p.Results.DefectLocations)
        % Generate random lattice positions (avoiding duplicates)
        defect_locations = zeros(num_defects, 2);
        placed_defects = 0;
        while placed_defects < num_defects
            % Generate random lattice coordinates
            new_loc = randi([1, N], 1, 2);
            
            % Check if this location is already used
            if placed_defects == 0 || ~any(all(defect_locations(1:placed_defects,:) == new_loc, 2))
                placed_defects = placed_defects + 1;
                defect_locations(placed_defects,:) = new_loc;
            end
        end
    else
        defect_locations = p.Results.DefectLocations;
        num_defects = size(defect_locations, 1);
    end
    
    % Create activation map (delta functions at defect locations)
    activation_map = zeros(grid_size, grid_size);
    grid_locations = round((defect_locations - 1) * scale_factor + 1);
    
    for i = 1:num_defects
        x = grid_locations(i,1);
        y = grid_locations(i,2);
        activation_map(y,x) = 1;
    end
    
    % Perform convolution
    multi_defect_LDoS = conv2(single_defect_LDoS, activation_map, 'same');
    
    % Visualization if requested
    if p.Results.PlotResult
        figure('Position', [100 100 1200 400]);
        
        % Single defect LDoS
        subplot(1,3,1);
        imagesc(single_defect_LDoS);
        title('Single Defect LDoS');
        colorbar;
        axis equal tight;
        xlabel('Grid X');
        ylabel('Grid Y');
        colormap(gca, gray);
        
        % Defect Locations Plot
        subplot(1,3,2);
        grid_img = zeros(grid_size, grid_size);
        for i = 1:size(grid_locations, 1)
            x = grid_locations(i,1);
            y = grid_locations(i,2);
            marker_size = 5;
            x_range = max(1,x-floor(marker_size/2)):min(grid_size,x+floor(marker_size/2));
            y_range = max(1,y-floor(marker_size/2)):min(grid_size,y+floor(marker_size/2));
            grid_img(y_range, x_range) = 1;
        end
        imagesc(grid_img);
        title(sprintf('%d Defect Locations', num_defects));
        colormap(gca, [1 1 1; 1 0 0]);  % White background, red markers
        axis equal tight;
        xlabel('Grid X');
        ylabel('Grid Y');
        
        % Add lattice grid lines
        hold on;
        for i = 1:N
            x_grid = (i-1) * scale_factor + 1;
            y_grid = (i-1) * scale_factor + 1;
            plot([1 grid_size], [y_grid y_grid], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'LineStyle', ':');
            plot([x_grid x_grid], [1 grid_size], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'LineStyle', ':');
        end
        hold off;
        
        % Convolved result
        subplot(1,3,3);
        imagesc(multi_defect_LDoS);
        title('Simulated Multi-defect LDoS');
        colorbar;
        axis equal tight;
        xlabel('Grid X');
        ylabel('Grid Y');
        colormap(gca, gray);
        
        sgtitle('Multi-defect LDoS Simulation', 'FontSize', 14);
    end
end