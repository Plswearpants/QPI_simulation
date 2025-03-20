function [horiz_data, diag_data] = extract_profiles(image_data, varargin)
% EXTRACT_PROFILES Extract horizontal and 45-degree profiles from 2D image
%
% Description:
%   Extracts both horizontal (0°) and diagonal (45°) profiles from an image
%   Returns raw and smoothed data for both profiles with position information
%
% Usage:
%   [horiz_data, diag_data] = extract_profiles(image_data)
%   [horiz_data, diag_data] = extract_profiles(image_data, options)
%
% Inputs:
%   image_data  - 2D matrix containing the image data
%   options     - Name-value pairs for additional parameters:
%       'center'      - [x,y] center point (default: center of image)
%       'width_pct'   - Percentage of full width to use (default: 80)
%       'bin_size'    - Bin size for directional mask (default: 2)
%       'bin_sep'     - Bin separation for directional mask (default: 1)
%       'units'       - 'pixels' or 'lattice' (default: 'pixels')
%       'grid_size'   - Grid size for lattice units (required if units='lattice')
%       'lattice_size'- Lattice size for conversion (required if units='lattice')
%
% Outputs:
%   horiz_data - Structure containing horizontal profile data:
%       .raw      - Raw intensity profile
%       .smooth   - Smoothed intensity profile
%       .x_raw    - Position values for raw profile
%       .x_smooth - Position values for smoothed profile
%   diag_data  - Structure containing diagonal (45°) profile data with the same fields

% Parse input arguments and pass them to theta_profile
p = inputParser;
addRequired(p, 'image_data', @isnumeric);
addParameter(p, 'center', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 2));
addParameter(p, 'width_pct', 80, @(x) isnumeric(x) && x > 0 && x <= 100);
addParameter(p, 'bin_size', 2, @(x) isnumeric(x) && x > 0);
addParameter(p, 'bin_sep', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'units', 'pixels', @(x) ischar(x) && any(strcmp(x, {'pixels', 'lattice'})));
addParameter(p, 'grid_size', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'lattice_size', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
parse(p, image_data, varargin{:});

% Extract horizontal (0 degrees) profile
[horiz_raw, horiz_smooth, horiz_x_raw, horiz_x_smooth] = theta_profile(image_data, 0, varargin{:});

% Extract diagonal (45 degrees) profile
[diag_raw, diag_smooth, diag_x_raw, diag_x_smooth] = theta_profile(image_data, 45, varargin{:});

% Package results in structures
horiz_data = struct('raw', horiz_raw, 'smooth', horiz_smooth, ...
                   'x_raw', horiz_x_raw, 'x_smooth', horiz_x_smooth);
diag_data = struct('raw', diag_raw, 'smooth', diag_smooth, ...
                  'x_raw', diag_x_raw, 'x_smooth', diag_x_smooth);

end 