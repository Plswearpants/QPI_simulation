%% load data and then select the slice to plot via d3gridDisplay
% Load the LDoS data
% Let user select the LDoS data file
[filename, pathname] = uigetfile({'*.mat', 'MAT-files (*.mat)'; '*.*', 'All Files (*.*)'}, 'Select LDoS data file');
if isequal(filename, 0) || isequal(pathname, 0)
    error('File selection canceled');
else
    fullpath = fullfile(pathname, filename);
    load(fullpath);
    fprintf('Loaded data from: %s\n', fullpath);
end
%% select the slice to process via d3gridDisplay
if ndims(LDoS_result) == 3
    % Display the 3D data detected message
    fprintf('3D LDoS data detected with %d energy slices\n', size(LDoS_result, 3));
    
    % Ask user for range type preference
    rangeType = questdlg('Select visualization range type:', 'Range Selection', ...
                         'global', 'dynamic', 'global');
    
    % First display the 3D data visualization
    figure; 
    d3gridDisplay_QPISIM(LDoS_result, rangeType);
    
    % Ask user if they want to select by slice number or energy
    selectionType = questdlg('Select slice by:', 'Selection Method', ...
                            'Slice Number', 'Energy Value', 'Slice Number');
    
    if strcmp(selectionType, 'Slice Number')
        % Let user select a slice using command line input
        fprintf('Enter slice number (1-%d) to process: ', size(LDoS_result, 3));
        sliceNum = input('');
        
        % Validate slice number
        if isempty(sliceNum) || ~isnumeric(sliceNum) || sliceNum < 1 || sliceNum > size(LDoS_result, 3)
            warning('Invalid slice number. Using slice 1.');
            sliceNum = 1;
        end
    else
        % Let user select by energy value
        fprintf('Available energy values (eV): \n');
        disp(omega_values);
        fprintf('Enter energy value to process: ');
        energyVal = input('');
        
        % Find closest energy value in omega_values
        [~, sliceNum] = min(abs(omega_values - energyVal));
        fprintf('Selected closest energy value: %.4f eV (slice %d)\n', omega_values(sliceNum), sliceNum);
    end
    
    % Extract the selected slice
    LDoS_noisy = LDoS_result(:,:,sliceNum);
    fprintf('Selected slice %d for processing (energy: %.4f eV)\n', sliceNum, omega_values(sliceNum));
else
    % If data is already 2D, use it directly
    LDoS_noisy = LDoS_result;
    fprintf('2D LDoS data detected\n');
    sliceNum = 1;
end

data = LDoS_result(:,:,sliceNum);
slice_energy = omega_values(sliceNum);

%% Vanila plot for the selected slice
figure;
imagesc(data);
colormap('gray');  % Use gray colormap for better contrast with green
axis square 
colorbar;
title(['\delta\rho @ ', num2str(slice_energy), ' eV']);

%% Plot the slice with lattice overlay
figure;
imagesc(data);
colormap('gray');  % Use gray colormap for better contrast with green
axis square;
colorbar;
title(['\delta\rho @ ', num2str(slice_energy), ' eV']);
hold on;

% Create lattice overlay
% N is already loaded from the data file
gridSize = size(data, 1);  % Get the size of the data grid
x_range = linspace(1, gridSize, N);  % Create N evenly spaced points in x
y_range = linspace(1, gridSize, N);  % Create N evenly spaced points in y

% Create the lattice grid
[X, Y] = meshgrid(x_range, y_range);

% Only keep points on the right half (including the middle line)
midpoint = 1;

right_half_mask = X >= midpoint;

% Plot the lattice points as green dots (only right half)
plot(X(right_half_mask), Y(right_half_mask), 'red.', 'MarkerSize', 5);

hold off;

%% Plot the horizontal and 45 degree profiles of the slice using maskDirectional
figure;

% Define default parameters for maskDirectional
profile_width = 3;  % Default width for the profile
bin_size = 2;       % Default bin size for smoothing
bin_sep = 1;        % Default bin separation for smoothing

% Create horizontal profile using maskDirectional
% Define start and end points for horizontal line through the middle
% Using 80% of the full width to avoid boundary effects
middle_row = round(size(data, 1)/2);
grid_width = size(data, 2);
margin = round(grid_width * 0.1);  % 10% margin on each side (80% length)
horiz_start = [margin + 1, middle_row];
horiz_end = [grid_width - margin, middle_row];

% Create masks for horizontal profile
[horiz_masks, horiz_masks_combined] = maskDirectional(data, ...
    'startPoint', horiz_start, ...
    'endPoint', horiz_end, ...
    'bin_size', bin_size, ...
    'bin_sep', bin_sep);

% Apply masks to extract horizontal profile
horiz_profile_raw = zeros(1, size(horiz_masks, 3));
for i = 1:size(horiz_masks, 3)
    mask = horiz_masks(:,:,i);
    horiz_profile_raw(i) = mean(data(mask));
end

% Apply combined masks for smoothed profile
horiz_profile_smooth = zeros(1, size(horiz_masks_combined, 3));
for i = 1:size(horiz_masks_combined, 3)
    mask = horiz_masks_combined(:,:,i);
    horiz_profile_smooth(i) = mean(data(mask));
end

% Create diagonal profile using maskDirectional
% Define start and end points for 45-degree diagonal (80% of full diagonal)
grid_min = min(size(data));
margin_diag = round(grid_min * 0.1);  % 10% margin on each side
diag_start = [margin_diag + 1, margin_diag + 1];
diag_end = [grid_min - margin_diag, grid_min - margin_diag];

% Create masks for diagonal profile
[diag_masks, diag_masks_combined] = maskDirectional(data, ...
    'startPoint', diag_start, ...
    'endPoint', diag_end, ...
    'bin_size', bin_size, ...
    'bin_sep', bin_sep);

% Apply masks to extract diagonal profile
diag_profile_raw = zeros(1, size(diag_masks, 3));
for i = 1:size(diag_masks, 3)
    mask = diag_masks(:,:,i);
    diag_profile_raw(i) = mean(data(mask));
end

% Apply combined masks for smoothed profile
diag_profile_smooth = zeros(1, size(diag_masks_combined, 3));
for i = 1:size(diag_masks_combined, 3)
    mask = diag_masks_combined(:,:,i);
    diag_profile_smooth(i) = mean(data(mask));
end

% Create subplot for horizontal profile
subplot(2, 1, 1);
% Convert pixel positions to lattice numbers using the ratio of grid size to lattice size
grid_to_lattice_ratio = gridSize / N; % Assuming gridSize and N are defined earlier
horiz_x = ((1:length(horiz_profile_raw)) - ceil(length(horiz_profile_raw)/2)) / grid_to_lattice_ratio;
plot(horiz_x, horiz_profile_raw, 'b-', 'LineWidth', 1, 'DisplayName', 'Raw');
hold on;
horiz_smooth_x = ((1:length(horiz_profile_smooth)) - ceil(length(horiz_profile_smooth)/2)) / grid_to_lattice_ratio;
plot(horiz_smooth_x, horiz_profile_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
title(['Horizontal Profile @ y = ', num2str(middle_row/grid_to_lattice_ratio), ' (80% width, avoiding boundaries)']);
xlabel('Lattice position');
ylabel('\delta\rho');
grid on;
legend('Location', 'best');

% Create subplot for diagonal profile
subplot(2, 1, 2);
% Convert pixel positions to lattice numbers using the ratio
diag_x = ((1:length(diag_profile_raw)) - ceil(length(diag_profile_raw)/2)) / grid_to_lattice_ratio;
plot(diag_x, diag_profile_raw, 'b-', 'LineWidth', 1, 'DisplayName', 'Raw');
hold on;
diag_smooth_x = ((1:length(diag_profile_smooth)) - ceil(length(diag_profile_smooth)/2)) / grid_to_lattice_ratio;
plot(diag_smooth_x, diag_profile_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
title(['45° Diagonal Profile (80% length, avoiding boundaries)']);
xlabel('Lattice position');
ylabel('\delta\rho');
grid on;
legend('Location', 'best');

sgtitle(['Profiles for \delta\rho @ ', num2str(slice_energy), ' eV', ...
    ' (bin\_size=', num2str(bin_size), ', bin\_sep=', num2str(bin_sep), ')']);

%% Start QPI process: 
% Generate QPI from target_LDoS and plot 
target_LDoS= LDoS_result;

QPI_sim= zeros(size(target_LDoS));
for k=1:size(target_LDoS,3)
    QPI_sim(:,:,k)=abs(fftshift(fft2(target_LDoS(:,:,k) - mean(mean(target_LDoS(:,:,k))))));
end
figure;
d3gridDisplay(QPI_sim,'dynamic')

