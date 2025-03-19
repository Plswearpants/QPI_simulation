% Script to compare multi-defect LDoS with convolution of single defect LDoS
% Author: Dong Chen
% Date: Feb 2025

%% Parameters and Data Loading
% Grid and lattice parameters

% Load your data here
% TODO: Replace these with your actual data loading commands
load('LDoS_single_defect.mat', 'LDoS_result','N');  % Load single defect LDoS
single_defect_LDoS = LDoS_result(:,:,1);
N_single = N;
single_grid_size = size(single_defect_LDoS,1);

load('LDoS_multi_defect', 'LDoS_result', 'used_locations', 'N');   % Load actual multi-defect LDoS
defect_locations = used_locations;
N_multi = N;
multi_defect_LDoS = LDoS_result(:,:,1);
grid_size = size(multi_defect_LDoS,1);


%% Compute Convolved LDoS
convolved_LDoS = compareDefectLDoS(single_defect_LDoS, defect_locations, N_multi, N_single, grid_size);

%% Visualization
figure('Position', [100 100 1000 800]);

% Single defect LDoS
subplot(2,2,1);
imagesc(single_defect_LDoS);
title('Single Defect LDoS');
colorbar;
axis equal tight;
xlabel('Grid X');
ylabel('Grid Y');
colormap(gca, gray); % Set to grayscale

% Defect Locations Plot
subplot(2,2,2);
% Create empty grid
grid_img = zeros(grid_size, grid_size);
% Convert lattice coordinates to grid coordinates
scale_factor = (grid_size - 1) / (N_multi - 1);
grid_locations = round((defect_locations - 1) * scale_factor + 1);
% Mark defect locations
for i = 1:size(grid_locations, 1)
    x = grid_locations(i,1);
    y = grid_locations(i,2);
    % Create a marker (5x5 square) at each defect location
    marker_size = 5;
    x_range = max(1,x-floor(marker_size/2)):min(grid_size,x+floor(marker_size/2));
    y_range = max(1,y-floor(marker_size/2)):min(grid_size,y+floor(marker_size/2));
    grid_img(y_range, x_range) = 1;
end
imagesc(grid_img);
title('Defect Locations');
colormap(gca, [1 1 1; 1 0 0]);  % White background, red markers
axis equal tight;
xlabel('Grid X');
ylabel('Grid Y');

% Add lattice grid lines with lighter color
hold on;
for i = 1:N_multi
    x_grid = (i-1) * scale_factor + 1;
    y_grid = (i-1) * scale_factor + 1;
    % Use a very light gray color for grid lines
    plot([1 grid_size], [y_grid y_grid], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'LineStyle', ':');
    plot([x_grid x_grid], [1 grid_size], 'Color', [0.8 0.8 0.8], 'LineWidth', 0.5, 'LineStyle', ':');
end
hold off;

% Convolved result
subplot(2,2,3);
imagesc(convolved_LDoS);
title('Convolved Multi-defect LDoS');
colorbar;
axis equal tight;
xlabel('Grid X');
ylabel('Grid Y');
colormap(gca, gray); % Set to grayscale

% Actual multi-defect LDoS
subplot(2,2,4);
imagesc(multi_defect_LDoS);
title('Actual Multi-defect LDoS');
colorbar;
axis equal tight;
xlabel('Grid X');
ylabel('Grid Y');
colormap(gca, gray); % Set to grayscale

% Add overall title
sgtitle('Multi-defect LDoS Comparison', 'FontSize', 14);

% Save results
%save('convolved_results.mat', 'convolved_LDoS');

% Print some statistics
correlation = corrcoef(multi_defect_LDoS(:), convolved_LDoS(:));
fprintf('Correlation coefficient between actual and convolved: %.4f\n', correlation(1,2));