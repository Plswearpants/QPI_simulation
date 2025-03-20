function [profile_raw, profile_smooth, x_raw, x_smooth] = theta_profile(image_data, theta_degrees, varargin)
% THETA_PROFILE Extract intensity profile along any angle in a 2D image
%
% Description:
%   Extracts intensity profiles along a specified angle from the center of the image
%   Can provide both raw and smoothed profiles with corresponding position values
%
% Usage:
%   [profile_raw, profile_smooth, x_raw, x_smooth] = theta_profile(image_data, theta_degrees)
%   [profile_raw, profile_smooth, x_raw, x_smooth] = theta_profile(image_data, theta_degrees, options)
%
% Inputs:
%   image_data     - 2D matrix containing the image data
%   theta_degrees  - Angle in degrees (0 = horizontal, 90 = vertical)
%   options        - Name-value pairs for additional parameters:
%       'center'      - [x,y] center point (default: center of image)
%       'width_pct'   - Percentage of full width to use (default: 80)
%       'bin_size'    - Bin size for directional mask (default: 2)
%       'bin_sep'     - Bin separation for directional mask (default: 1)
%       'units'       - 'pixels' or 'lattice' (default: 'pixels')
%       'grid_size'   - Grid size for lattice units (required if units='lattice')
%       'lattice_size'- Lattice size for conversion (required if units='lattice')
%
% Outputs:
%   profile_raw    - Raw intensity profile along the specified angle
%   profile_smooth - Smoothed intensity profile along the specified angle
%   x_raw          - Position values for raw profile (centered at 0)
%   x_smooth       - Position values for smoothed profile (centered at 0)

% Parse input arguments
p = inputParser;
addRequired(p, 'image_data', @isnumeric);
addRequired(p, 'theta_degrees', @isnumeric);
addParameter(p, 'center', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 2));
addParameter(p, 'width_pct', 80, @(x) isnumeric(x) && x > 0 && x <= 100);
addParameter(p, 'bin_size', 2, @(x) isnumeric(x) && x > 0);
addParameter(p, 'bin_sep', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'units', 'pixels', @(x) ischar(x) && any(strcmp(x, {'pixels', 'lattice'})));
addParameter(p, 'grid_size', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'lattice_size', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
parse(p, image_data, theta_degrees, varargin{:});

center = p.Results.center;
width_pct = p.Results.width_pct;
bin_size = p.Results.bin_size;
bin_sep = p.Results.bin_sep;
units = p.Results.units;
grid_size = p.Results.grid_size;
lattice_size = p.Results.lattice_size;

% Check if lattice conversion is requested but parameters are missing
if strcmp(units, 'lattice') && (isempty(grid_size) || isempty(lattice_size))
    error('For lattice units, both grid_size and lattice_size must be provided');
end

% Get image dimensions
[img_height, img_width] = size(image_data);

% Set default center if not provided
if isempty(center)
    center = [ceil(img_width/2), ceil(img_height/2)];
end

% Convert angle to radians
theta_rad = theta_degrees * pi / 180;

% Calculate the distance from center to edge along the angle
if theta_degrees == 0 || theta_degrees == 180
    % Horizontal line
    dist_to_edge = max(center(1) - 1, img_width - center(1));
elseif theta_degrees == 90 || theta_degrees == 270
    % Vertical line
    dist_to_edge = max(center(2) - 1, img_height - center(2));
else
    % Angled line - calculate distance to edges
    dx = cos(theta_rad);
    dy = sin(theta_rad);
    
    % Calculate distances to edges
    if dx > 0
        tx1 = (img_width - center(1)) / dx;
    elseif dx < 0
        tx1 = (1 - center(1)) / dx;
    else
        tx1 = Inf;
    end
    
    if dy > 0
        ty1 = (img_height - center(2)) / dy;
    elseif dy < 0
        ty1 = (1 - center(2)) / dy;
    else
        ty1 = Inf;
    end
    
    % Use minimum distance to determine edge point
    dist_to_edge = min(tx1, ty1);
end

% Apply width percentage to limit the distance
dist_to_use = dist_to_edge * (width_pct / 100);

% Calculate start and end points
start_x = round(center(1) - dist_to_use * cos(theta_rad));
start_y = round(center(2) - dist_to_use * sin(theta_rad));
end_x = round(center(1) + dist_to_use * cos(theta_rad));
end_y = round(center(2) + dist_to_use * sin(theta_rad));

% Ensure coordinates are within image bounds
start_x = max(1, min(img_width, start_x));
start_y = max(1, min(img_height, start_y));
end_x = max(1, min(img_width, end_x));
end_y = max(1, min(img_height, end_y));

% Create masks for profile
try
    [masks, masks_combined] = maskDirectional(image_data, ...
        'startPoint', [start_x, start_y], ...
        'endPoint', [end_x, end_y], ...
        'bin_size', bin_size, ...
        'bin_sep', bin_sep);
    
    % Apply masks to extract profile
    profile_raw = zeros(1, size(masks, 3));
    for i = 1:size(masks, 3)
        mask = masks(:,:,i);
        profile_raw(i) = mean(image_data(mask));
    end
    
    % Apply combined masks for smoothed profile
    profile_smooth = zeros(1, size(masks_combined, 3));
    for i = 1:size(masks_combined, 3)
        mask = masks_combined(:,:,i);
        profile_smooth(i) = mean(image_data(mask));
    end
    
catch e
    warning('Error creating profile with maskDirectional: %s\nUsing direct extraction instead.', e.message);
    
    % Fallback to direct line extraction using Bresenham's algorithm
    [x_line, y_line] = bresenham_line(start_x, start_y, end_x, end_y);
    
    % Extract values along the line
    valid_indices = x_line >= 1 & x_line <= img_width & y_line >= 1 & y_line <= img_height;
    x_line = x_line(valid_indices);
    y_line = y_line(valid_indices);
    
    % Convert indices to linear indices
    lin_indices = sub2ind(size(image_data), y_line, x_line);
    profile_raw = image_data(lin_indices);
    
    % Simple smoothing for the fallback method
    profile_smooth = smoothdata(profile_raw, 'gaussian', min(5, floor(length(profile_raw)/4)));
end

% Calculate position values
if strcmp(units, 'pixels')
    % Use pixels centered at 0
    x_raw = linspace(-length(profile_raw)/2, length(profile_raw)/2, length(profile_raw));
    x_smooth = linspace(-length(profile_smooth)/2, length(profile_smooth)/2, length(profile_smooth));
else
    % Convert to lattice units
    grid_to_lattice_ratio = grid_size / lattice_size;
    x_raw = ((1:length(profile_raw)) - ceil(length(profile_raw)/2)) / grid_to_lattice_ratio;
    x_smooth = ((1:length(profile_smooth)) - ceil(length(profile_smooth)/2)) / grid_to_lattice_ratio;
    
    % If this is a diagonal, apply the appropriate scaling factor
    if mod(theta_degrees, 90) ~= 0
        % For non-cardinal angles, scale by the diagonal factor
        diagonal_factor = sqrt(1 + tan(theta_rad)^2);
        x_raw = x_raw * diagonal_factor;
        x_smooth = x_smooth * diagonal_factor;
    end
end

end

% Helper function for Bresenham's line algorithm
function [x, y] = bresenham_line(x1, y1, x2, y2)
    % Implementation of Bresenham's line algorithm to get all points along a line
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    steep = dy > dx;
    
    if steep
        % Swap x and y
        [x1, y1] = deal(y1, x1);
        [x2, y2] = deal(y2, x2);
    end
    
    if x1 > x2
        % Swap start and end
        [x1, x2] = deal(x2, x1);
        [y1, y2] = deal(y2, y1);
    end
    
    dx = x2 - x1;
    dy = abs(y2 - y1);
    error = dx / 2;
    
    if y1 < y2
        ystep = 1;
    else
        ystep = -1;
    end
    
    y = y1;
    x_points = zeros(1, x2-x1+1);
    y_points = zeros(1, x2-x1+1);
    
    for i = 1:(x2-x1+1)
        x = x1 + i - 1;
        if steep
            x_points(i) = y;
            y_points(i) = x;
        else
            x_points(i) = x;
            y_points(i) = y;
        end
        
        error = error - dy;
        if error < 0
            y = y + ystep;
            error = error + dx;
        end
    end
    
    x = x_points;
    y = y_points;
end 