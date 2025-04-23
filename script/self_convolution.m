% Self-convolution demonstration script
% This script demonstrates self-convolution of hollow shapes (square and circle)
% and visualizes both the original images and their self-convolutions

%% Parameters
gridSize = 500;
square_size = 150;
thickness = 5;  % 2-pixel thickness
load('InverseGray', 'invgray');
%% First figure: Original Hollow Square
figure('Name', 'Original Hollow Square', 'Position', [100, 100, 400, 400]);
square = zeros(gridSize, gridSize);
start_idx = round((gridSize - square_size)/2);
end_idx = start_idx + square_size - 1;

% Create outer square
square(start_idx:end_idx, start_idx:end_idx) = 1;
% Remove inner square
inner_start = start_idx + thickness;
inner_end = end_idx - thickness;
square(inner_start:inner_end, inner_start:inner_end) = 0;

imagesc(square);
axis equal tight;
colormap(invgray);
title('Original Hollow Square');
xlabel('x (pixels)');
ylabel('y (pixels)');
colorbar;

%% Prepare Circle for Convolution
% Calculate perimeter of square for circle radius
square_perimeter = 4 * square_size;
radius = square_perimeter/(2*pi);
circle = zeros(gridSize, gridSize);
center = round(gridSize/2);
[X, Y] = meshgrid(1:gridSize, 1:gridSize);
distances = sqrt((X - center).^2 + (Y - center).^2);
circle = (distances <= radius) & (distances > (radius - thickness));

%% Compute both convolutions and determine common scale
% Create circular mask for center region
mask = sqrt((X - center).^2 + (Y - center).^2) > thickness;

% Compute convolutions
square_conv = conv2(square, square, 'same');
circle_conv = conv2(circle, circle, 'same');

% Convert both to log scale
square_conv_log = log10(square_conv + 1);
circle_conv_log = log10(circle_conv + 1);

% Find common min and max for normalization
min_val = min(min(square_conv_log(:)), min(circle_conv_log(:)));
max_val = max(max(square_conv_log(:)), max(circle_conv_log(:)));

% Normalize both using the same scale
square_conv_norm = (square_conv_log - min_val) / (max_val - min_val);
circle_conv_norm = (circle_conv_log - min_val) / (max_val - min_val);

%% Second figure: Square Convolution
figure('Name', 'Square Convolution', 'Position', [550, 100, 400, 400]);
imagesc(square_conv_norm);
axis equal tight;
colormap(invgray);
title('Log of Self-Convolution (Square)');
xlabel('x (pixels)');
ylabel('y (pixels)');
c = colorbar;
%c.Label.String = 'Normalized Intensity (log scale)';
caxis([0 1]);  % Set same color axis limits

%% Third figure: Original Hollow Circle
figure('Name', 'Original Hollow Circle', 'Position', [100, 550, 400, 400]);
imagesc(circle);
axis equal tight;
colormap(invgray);
title(sprintf('Original Hollow Circle (r = %.2f)', radius));
xlabel('x (pixels)');
ylabel('y (pixels)');
colorbar;

%% Fourth figure: Circle Convolution
figure('Name', 'Circle Convolution', 'Position', [550, 550, 400, 400]);
imagesc(circle_conv_norm);
axis equal tight;
colormap(invgray);
title('Log of Self-Convolution (Circle)');
xlabel('x (pixels)');
ylabel('y (pixels)');
colorbar;

