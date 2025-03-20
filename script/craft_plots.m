%% Load data
[filename, pathname] = uigetfile({'*.mat'}, 'Select LDoS data file');
if isequal(filename, 0)
    error('File selection canceled');
end
load(fullfile(pathname, filename));
fprintf('Loaded data from: %s\n', fullfile(pathname, filename));

%% Load inversegray colormap to use consistently throughout the script
if ~exist('invgray', 'var')
    try
        load('InverseGray', 'invgray');
    catch
        % Create a simple inverse grayscale colormap if file not found
        invgray = 1 - gray(256);
        warning('InverseGray.mat not found. Created a simple inverse grayscale colormap.');
    end
end

%% Get user preferences
config = get_user_preferences();

%% Process data
data = process_ldos_data(LDoS_result, omega_values, config);
gridSize = size(data.slice, 1);

%% Plot LDoS slice with lattice overlay
plot_slice_with_lattice(data.slice, N, gridSize, config.axisMode, data.energy);

%% Plot LDoS profiles
figure;
[horiz_data, diag_data] = extract_profiles(data.slice, ...
    'bin_size', config.bin_size, ...
    'bin_sep', config.bin_sep, ...
    'width_pct', config.width_pct, ...
    'units', 'lattice', ...
    'grid_size', gridSize, ...
    'lattice_size', N);

% Plot horizontal profile
subplot(2, 1, 1);
plot(horiz_data.x_raw, horiz_data.raw, 'b-', 'LineWidth', 1, 'DisplayName', 'Raw');
hold on;
plot(horiz_data.x_smooth, horiz_data.smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
title(['Horizontal Profile (', num2str(config.width_pct), '% width)']);
xlabel('Lattice position'); ylabel('\delta\rho');
grid on; legend('Location', 'best');

% Plot diagonal profile
subplot(2, 1, 2);
plot(diag_data.x_raw, diag_data.raw, 'b-', 'LineWidth', 1, 'DisplayName', 'Raw');
hold on;
plot(diag_data.x_smooth, diag_data.smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
title(['45° Diagonal Profile (', num2str(config.width_pct), '% length)']);
xlabel('Lattice position'); ylabel('\delta\rho');
grid on; legend('Location', 'best');

sgtitle(['Profiles for \delta\rho @ ', num2str(data.energy), ' eV']);

%% Generate and display QPI
QPI_sim = generate_qpi(data.used);
figure;
d3gridDisplay(QPI_sim, 'dynamic');
title('Full QPI');

% Ask users if they want to crop the QPI signal
crop_qpi = questdlg('Do you want to crop the QPI signal to focus on the central region between Bragg peaks?', 'Crop QPI', ...
    'Yes', 'No', 'Yes');

% Get slice indices for QPI visualization
[slice_indices, slice_energies] = get_qpi_slices(data, config);

% Compute the bragg peak positions
pixels_per_unit_cell = gridSize / N;
edge2bragg_peak = 2*round(gridSize/ pixels_per_unit_cell);

% Initialize QPI data for plotting
QPI_for_plotting = QPI_sim;

% Crop if user selected yes
if strcmp(crop_qpi, 'Yes')
    % Crop the QPI signal to the central region between Bragg peaks
    QPI_sim_cropped = QPI_sim(edge2bragg_peak:end-edge2bragg_peak+1, edge2bragg_peak:end-edge2bragg_peak+1, :);
    
    % Display the cropped QPI
    figure('Name', 'Cropped QPI');
    d3gridDisplay(QPI_sim_cropped, 'dynamic');
    title('Cropped QPI (Central Region)');
    
    % Use cropped data for slice plotting
    QPI_for_plotting = QPI_sim_cropped;
end

% Plot QPI slices
plot_qpi_grid(QPI_for_plotting, slice_indices, data.omega_values, config.qpiAxisMode, N, gridSize);

%% Plot QPI profiles
figure;
for i = 1:length(slice_indices)
    slice_idx = slice_indices(i);
    qpi_slice = QPI_sim(:,:,slice_idx);
    
    % Extract profiles for QPI
    if strcmp(config.qpiAxisMode, 'Grid Pixels')
        [qpi_horiz, qpi_diag] = extract_profiles(qpi_slice, ...
            'bin_size', config.bin_size, ...
            'bin_sep', config.bin_sep, ...
            'width_pct', config.width_pct, ...
            'units', 'pixels');
    else
        [qpi_horiz, qpi_diag] = extract_profiles(qpi_slice, ...
            'bin_size', config.bin_size, ...
            'bin_sep', config.bin_sep, ...
            'width_pct', config.width_pct, ...
            'units', 'lattice', ...
            'grid_size', gridSize, ...
            'lattice_size', N);
    end
    
    % Plot profiles (horizontal and diagonal)
    subplot(2, length(slice_indices), i);
    plot_profile(qpi_horiz, config.qpiAxisMode, slice_energies(i), 'Horizontal');
    
    subplot(2, length(slice_indices), i + length(slice_indices));
    plot_profile(qpi_diag, config.qpiAxisMode, slice_energies(i), 'Diagonal');
end

sgtitle(['QPI Intensity Profiles (bin\_size=', num2str(config.bin_size), ...
    ', bin\_sep=', num2str(config.bin_sep), ')']);

%% Create and plot band structure

% Create the band structure for CEC calculation 
if ~exist('E', 'var') || ~exist('kx', 'var') || ~exist('ky', 'var')
    % Define parameters for the dispersion relation if not already defined
    a = 1e-9; % lattice constant
    t = -0.2; % hopping parameter
    E0 = 0; % on-site energy
    
    % Define k-space grid
    k_vals = linspace(-pi/a, pi/a, n);
    [kx, ky] = meshgrid(k_vals, k_vals);
    
    % Compute the energy dispersion
    E = E0 - 2 * t * (cos(kx * a) + cos(ky * a));
end

% Plot the energy dispersion
figure;
surf(kx, ky, E);
colormap("jet");
colorbar;
title('Energy Dispersion Relation');
xlabel('k_x');
ylabel('k_y');
zlabel('Energy (E)');
grid on;

%% Create CEC and JDOS(self-conv of CEC) with finer energy resolution
omega_vals = linspace(min(E(:)), max(E(:)), 41); % Keep 41 slices for visualization
energy_tolerance = 0.01; % Much finer energy broadening (0.001 eV)

% Initialize arrays to store CEC and convolution results
CEC = zeros(size(kx,1), size(kx,2), length(omega_vals));
CEC_conv = zeros(2*size(kx,1), 2*size(kx,2), length(omega_vals));

% Create expanded k-space grid for convolution (double the range)
k_vals_expanded = linspace(-2*pi/a, 2*pi/a, 2*size(kx,1));
[kx_expanded, ky_expanded] = meshgrid(k_vals_expanded, k_vals_expanded);

% Compute CEC and convolution for each energy
for i = 1:length(omega_vals)
    % Create constant energy contour with much smaller tolerance (0.001 eV)
    CEC(:,:,i) = double(abs(E - omega_vals(i)) < energy_tolerance);
    
    % Check if CEC is too sparse (can happen with very small tolerance)
    if sum(CEC(:,:,i), 'all') < 10
        warning(['Very few points in CEC at energy ' num2str(omega_vals(i)) ...
                 '. Consider increasing k-space resolution or energy tolerance.']);
        
        % Option: Adaptively increase tolerance if needed
        adaptive_tolerance = energy_tolerance;
        while sum(CEC(:,:,i), 'all') < 10 && adaptive_tolerance < 0.01
            adaptive_tolerance = adaptive_tolerance * 2;
            CEC(:,:,i) = double(abs(E - omega_vals(i)) < adaptive_tolerance);
        end
        
        fprintf('Used adaptive tolerance of %.6f eV for energy %.4f eV\n', ...
                adaptive_tolerance, omega_vals(i));
    end
    
    % Compute convolution properly with FFT method to preserve correct k-space scaling
    fft_CEC = fft2(CEC(:,:,i), 2*size(CEC,1), 2*size(CEC,2)); % Zero-pad for expanded k-space
    conv_result = abs(ifft2(fft_CEC .* fft_CEC)); % Convolution in real space is multiplication in Fourier space
    CEC_conv(:,:,i) = conv_result / max(conv_result(:), [], 'all'); % Normalize
end

%% Visualize CEC (Constant Energy Contours)
figure('Name', 'Constant Energy Contours (CEC)');
d3gridDisplay(CEC, 'dynamic');
title('Constant Energy Contours (CEC)');

%% Visualize JDOS (Joint Density of States)
figure('Name', 'Joint Density of States (JDOS)');
d3gridDisplay(CEC_conv, 'dynamic');
title('Joint Density of States (JDOS)');

%% CEC and JDOS slice visualization with reciprocal space units
% Get the single chosen slice index
chosen_slice = slice_indices(1); % Take only the first index if multiple were selected

% Set up the reciprocal space grid
k_vals = linspace(-pi/a, pi/a, size(kx, 1)); % Original k-space grid for CEC
k_vals_norm = k_vals / (pi/a); % Normalize to units of π/a

% For JDOS (which has double the k-space range)
k_vals_expanded = linspace(-2*pi/a, 2*pi/a, size(CEC_conv, 1));
k_vals_expanded_norm = k_vals_expanded / (pi/a); % Normalize to units of π/a

% Apply dynamic range to CEC
figure('Name', ['CEC at E = ' num2str(omega_values(chosen_slice)) ' eV']);
cec_slice = CEC(:,:,chosen_slice);

% Calculate dynamic range parameters (median ± nos*std)
nos = 25; % Controls contrast
med_val = median(cec_slice(:));
std_val = std(cec_slice(:));
dyn_range = [med_val-nos*std_val, med_val+nos*std_val];

% Display with invgray colormap, dynamic range, and proper k-space axes
imagesc(k_vals_norm, k_vals_norm, cec_slice, dyn_range);
colormap(invgray);
colorbar;
axis square;
xlabel('k_x (\pi/a)');
ylabel('k_y (\pi/a)');
title(['CEC at E = ' num2str(omega_values(chosen_slice)) ' eV ']);

% Apply dynamic range to JDOS
figure('Name', ['JDOS at E = ' num2str(omega_values(chosen_slice)) ' eV']);
jdos_slice = CEC_conv(:,:,chosen_slice);

% Calculate dynamic range parameters for JDOS
med_val_jdos = median(jdos_slice(:));
std_val_jdos = std(jdos_slice(:));
dyn_range_jdos = [med_val_jdos-nos*std_val_jdos, med_val_jdos+nos*std_val_jdos];

% Display with invgray colormap, dynamic range, and proper k-space axes
imagesc(k_vals_expanded_norm, k_vals_expanded_norm, jdos_slice, dyn_range_jdos);
colormap(invgray);
colorbar;
axis square;
xlabel('k_x (\pi/a)');
ylabel('k_y (\pi/a)');
title(['JDOS at E = ' num2str(omega_values(chosen_slice)) ' eV']);

%% QPI cropped visualization with same slice and dynamic range
% Use the same chosen slice as above
figure('Name', ['QPI (cropped) at E = ' num2str(omega_values(chosen_slice)) ' eV']);

% Get the QPI data for the chosen slice
qpi_slice = QPI_for_plotting(:,:,chosen_slice);

% Use the same dynamic range approach as for CEC and JDOS
med_val_qpi = median(qpi_slice(:));
std_val_qpi = std(qpi_slice(:));
dyn_range_qpi = [med_val_qpi-nos*std_val_qpi, med_val_qpi+nos*std_val_qpi];

% Display with proper k-space axes using the same normalization as JDOS
pixels_per_unit_cell = gridSize / N;
zoom_ratio = (size(QPI_for_plotting,1)/gridSize);
kmax = pi * pixels_per_unit_cell/2 * zoom_ratio;
k_range = linspace(-kmax, kmax, size(QPI_for_plotting,1)) / pi;

% Display the QPI with the same colormap and dynamic range
imagesc(k_range, k_range, qpi_slice, dyn_range_qpi);
colormap(invgray);
colorbar;
axis square;
xlabel('k_x (\pi/a)');
ylabel('k_y (\pi/a)');
title(['QPI at E = ' num2str(omega_values(chosen_slice)) ' eV']);


%% Helper functions 
function config = get_user_preferences()
    % Get all user preferences and settings at once
    config = struct();
    
    % Visualization preferences
    config.rangeType = questdlg('Select visualization range type:', 'Range Selection', ...
                             'global', 'dynamic', 'global');
    
    config.axisMode = questdlg('Select axis display mode:', 'Axis Mode', ...
                        'Grid Pixels', 'Physical Units', 'Grid Pixels');
    
    config.qpiAxisMode = questdlg('Select QPI axis display mode:', 'QPI Axis Mode', ...
                       'Grid Pixels', 'Reciprocal Space Units', 'Grid Pixels');
    
    % Profile settings
    config.bin_size = 2;      % Bin size for smoothing
    config.bin_sep = 1;       % Bin separation for smoothing
    config.width_pct = 80;    % Width percentage for profiles
    
    % Energy selection preferences
    config.input_type = questdlg('Select slices by:', 'Slice Selection Method', ...
                              'Energy Values', 'Slice Indices', 'Energy Values');
                          
    if strcmp(config.input_type, 'Energy Values')
        config.input_type = 'energy';
    else
        config.input_type = 'slice';
    end
    
    % Default values
    if strcmp(config.input_type, 'energy')
        config.slice_energies = [-0.05]; % Default energy value in eV
    else
        config.slice_indices_direct = [1, 5, 10, 15, 20, 25]; % Default slice indices
    end
end

function data = process_ldos_data(LDoS_result, omega_values, config)
    % Process the LDoS data, including masking defects if requested
    
    % Initialize data structure with original unprocessed data
    data = struct();
    data.original = LDoS_result;  % Store original unprocessed data
    
    % Handle defect masking
    if questdlg('Do you want to mask defects in the LDoS data?', 'Mask Defects', ...
               'Yes', 'No', 'Yes') == "Yes"
        
        centerSliceIdx = round(size(LDoS_result, 3) / 2);
        centerEnergy = omega_values(centerSliceIdx);
        
        fprintf('Please click on defects in the image at energy %.4f eV (center slice).\n', centerEnergy);
        
        if questdlg('Use predefined defect locations?', 'Defect Locations', ...
                   'Yes', 'No (Interactive)', 'No (Interactive)') == "Yes"
            defectLocations = input('Enter defect coordinates as [x1,y1; x2,y2; ...]: ');
            if ~isnumeric(defectLocations) || size(defectLocations, 2) ~= 2
                warning('Invalid input format. Switching to interactive mode.');
                defectLocations = [];
                LDoS_processed = defectMask(LDoS_result, centerEnergy, omega_values);
            else
                LDoS_processed = defectMask(LDoS_result, centerEnergy, omega_values, 'locations', defectLocations);
            end
        else
            LDoS_processed = defectMask(LDoS_result, centerEnergy, omega_values);
        end
        
        data.processed = LDoS_processed;  % Store processed data
        fprintf('Defect masking complete.\n');
    else
        % No processing was requested - processed data same as original
        data.processed = [];  % No separate processed data
        fprintf('Proceeding with original (unprocessed) data.\n');
    end
    
    % Choose the appropriate data for further operations
    % If processed data exists, use it; otherwise use original
    working_data = data.processed;
    if isempty(working_data)
        working_data = data.original;
    end
    
    % Store the working data for convenience
    data.used = working_data;
    
    % Select slice(s)
    if ndims(working_data) == 3
        fprintf('3D LDoS data detected with %d energy slices\n', size(working_data, 3));
        
        % Display data for browsing
        figure; 
        d3gridDisplay_QPISIM(working_data, config.rangeType);
        
        % Get list of slices to analyze
        if questdlg('Select multiple slices?', 'Multiple Slices', ...
                   'Yes', 'No (single slice)', 'No (single slice)') == "Yes"
            
            % Multiple slice selection mode
            if strcmp(config.input_type, 'energy')
                % Energy-based selection
                fprintf('Available energy values (eV):\n');
                disp(omega_values);
                
                % Get list of energies
                energy_list = input('Enter list of energy values [e1, e2, ...]: ');
                if ~isnumeric(energy_list) || isempty(energy_list)
                    warning('Invalid input. Using default energy value.');
                    energy_list = config.slice_energies;
                end
                
                % Find nearest indices for each energy
                slice_indices = zeros(size(energy_list));
                for i = 1:length(energy_list)
                    [~, slice_indices(i)] = min(abs(omega_values - energy_list(i)));
                    fprintf('Selected energy %.4f eV (slice %d)\n', energy_list(i), slice_indices(i));
                end
                
                % Store for later use
                config.slice_energies = energy_list;
                data.multi_energies = energy_list;
                data.multi_slices = slice_indices;
                
                % Use the first slice as the primary one for analysis
                sliceNum = slice_indices(1);
                
            else
                % Index-based selection
                fprintf('Enter slice indices (1-%d):\n', size(working_data, 3));
                
                % Get list of indices
                slice_list = input('Enter list of slice indices [i1, i2, ...]: ');
                if ~isnumeric(slice_list) || isempty(slice_list) || ...
                   any(slice_list < 1) || any(slice_list > size(working_data, 3))
                    warning('Invalid input. Using default slice indices.');
                    slice_list = min(config.slice_indices_direct, size(working_data, 3));
                end
                
                % Store for later use
                config.slice_indices_direct = slice_list;
                data.multi_slices = slice_list;
                data.multi_energies = omega_values(slice_list);
                
                % Use the first slice as the primary one for analysis
                sliceNum = slice_list(1);
            end
            
            fprintf('Multiple slices selected. Using slice %d (%.4f eV) as primary.\n', ...
                    sliceNum, omega_values(sliceNum));
            
        else
            % Single slice selection mode
            if questdlg('Select slice by:', 'Selection Method', ...
                   'Slice Number', 'Energy Value', 'Slice Number') == "Slice Number"
                sliceNum = input(['Enter slice number (1-', num2str(size(working_data, 3)), '): ']);
                
                if isempty(sliceNum) || ~isnumeric(sliceNum) || sliceNum < 1 || sliceNum > size(working_data, 3)
                    warning('Invalid slice number. Using slice 1.');
                    sliceNum = 1;
                end
            else
                disp('Available energy values (eV):');
                disp(omega_values);
                energyVal = input('Enter energy value: ');
                [~, sliceNum] = min(abs(omega_values - energyVal));
                fprintf('Selected closest energy value: %.4f eV (slice %d)\n', omega_values(sliceNum), sliceNum);
            end
            
            % Store single slice as multi-slice list for consistency
            data.multi_slices = sliceNum;
            data.multi_energies = omega_values(sliceNum);
        end
        
        % Store selected slice and related info
        data.slice = working_data(:,:,sliceNum);
        data.energy = omega_values(sliceNum);
    else
        % Handle 2D data
        data.slice = working_data;
        data.energy = NaN;
        sliceNum = 1;
        data.multi_slices = 1;
        data.multi_energies = NaN;
    end
    
    % Store additional metadata
    data.sliceNum = sliceNum;
    data.omega_values = omega_values;
    
    % Add a flag to indicate if processing was applied
    data.is_processed = ~isempty(data.processed);
end

function QPI = generate_qpi(LDoS_data)
    % Generate QPI from LDoS data
    QPI = zeros(size(LDoS_data));
    for k = 1:size(LDoS_data, 3)
        QPI(:,:,k) = abs(fftshift(fft2(LDoS_data(:,:,k) - mean(mean(LDoS_data(:,:,k))))));
    end
end

function plot_slice_with_lattice(data, N, gridSize, axisMode, energy)
    % Plot LDoS slice with lattice overlay
    figure;
    
    % Load invgray colormap from base workspace
    if ~exist('invgray', 'var')
        try
            invgray = evalin('base', 'invgray');
        catch
            % Create a simple inverse grayscale colormap if not found
            invgray = 1 - gray(256);
        end
    end
    
    if strcmp(axisMode, 'Grid Pixels')
        imagesc(data);
        xlabel('Pixel X');
        ylabel('Pixel Y');
        title(['\delta\rho @ ', num2str(energy), ' eV (Grid Pixels)']);
    else
        imagesc([1, N], [1, N], data);
        xlabel('Lattice Position X');
        ylabel('Lattice Position Y');
        title(['\delta\rho @ ', num2str(energy), ' eV (Lattice Units)']);
    end
    
    colormap(invgray); % Use inversegray colormap
    axis square;
    colorbar;
    
    hold on;
    
    % Create lattice overlay
    x_range = linspace(1, gridSize, N);
    y_range = linspace(1, gridSize, N);
    [X, Y] = meshgrid(x_range, y_range);
    
    right_half_mask = X >= 1;
    
    if strcmp(axisMode, 'Grid Pixels')
        plot(X(right_half_mask), Y(right_half_mask), 'red.', 'MarkerSize', 5);
    else
        lattice_X = linspace(1, N, length(x_range));
        lattice_Y = linspace(1, N, length(y_range));
        [Lattice_X, Lattice_Y] = meshgrid(lattice_X, lattice_Y);
        plot(Lattice_X(right_half_mask), Lattice_Y(right_half_mask), 'red.', 'MarkerSize', 5);
    end
    
    hold off;
end

function plot_qpi_grid(QPI, slice_indices, slice_energies, qpiAxisMode, N, gridSize)
    % Plot QPI slices in a grid layout
    n = length(slice_indices);
    
    % Determine layout
    if mod(n, 2) == 0
        rows = 2;
        cols = n/2;
    else
        rows = 1;
        cols = n;
    end
    
    figure('Position', [100, 100, 200*cols, 200*rows]);
    
    % Load invgray colormap from base workspace
    if ~exist('invgray', 'var')
        try
            invgray = evalin('base', 'invgray');
        catch
            % Create a simple inverse grayscale colormap if not found
            invgray = 1 - gray(256);
        end
    end
    
    for i = 1:n
        subplot(rows, cols, i);
        slice_idx = slice_indices(i);
        
        if strcmp(qpiAxisMode, 'Grid Pixels')
            imagesc(QPI(:,:,slice_idx));
            xlabel('Pixel X');
            ylabel('Pixel Y');
        else
            pixels_per_unit_cell = gridSize / N;
            zoom_ratio = (size(QPI,1)/gridSize);
            kmax = pi * pixels_per_unit_cell * zoom_ratio;
            k_range = linspace(-kmax, kmax, size(QPI,1)) / pi;
            
            imagesc(k_range, k_range, QPI(:,:,slice_idx));
            xlabel('k_x (pi/a)');
            ylabel('k_y (pi/a)');
            
            hold on;
            bragg_positions = [2, 0; 0, -2; -2, 0; 0, 2];
            scatter(bragg_positions(:,1), bragg_positions(:,2), 15, 'r', 'filled', 'o', 'MarkerEdgeColor', 'k');
            hold off;
        end
        
        colormap(invgray); % Use inversegray colormap
        axis square;
        colorbar;
        
        title(['E = ', num2str(slice_energies(slice_idx), '%.3f'), ' eV']);
    end
    
    if strcmp(qpiAxisMode, 'Grid Pixels')
        sgtitle('QPI at Selected Energies (Grid Pixels)');
    else
        sgtitle('QPI at Selected Energies (Reciprocal Space Units)');
    end
end

function plot_profile(profile_data, mode, energy, direction)
    % Helper function to plot profile data
    plot(profile_data.x_raw, profile_data.raw, 'b-', 'LineWidth', 1, 'DisplayName', 'Raw');
    hold on;
    plot(profile_data.x_smooth, profile_data.smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed');
    
    if strcmp(mode, 'Grid Pixels')
        xlabel('Distance from Center (pixels)');
    elseif strcmp(direction, 'Horizontal')
        xlabel('k_x (2\pi/a)');
    else
        xlabel('k_{diagonal} (2\pi/a)');
    end
    
    ylabel('Intensity');
    title([direction, ' Profile, E = ', num2str(energy, '%.2f'), ' eV']);
    grid on;
    legend('Location', 'best');
end

function [slice_indices, slice_energies] = get_qpi_slices(data, config)
    % Function to get QPI slice indices either from config or user input
    
    % Ask if user wants to use the same slices as LDoS or specify new ones
    use_same = questdlg('Use the same slices as LDoS visualization?', 'QPI Slices', ...
                        'Yes', 'No (specify new)', 'Yes');
                        
    if strcmp(use_same, 'Yes')
        % Use the same slices as LDoS
        slice_indices = data.multi_slices;
        slice_energies = data.multi_energies;
        fprintf('Using the same %d slice(s) as LDoS visualization.\n', length(slice_indices));
    else
        % Specify new slices for QPI
        if strcmp(config.input_type, 'energy')
            % Energy-based selection
            fprintf('Available energy values (eV):\n');
            disp(data.omega_values);
            
            % Get list of energies
            energy_list = input('Enter list of energy values for QPI [e1, e2, ...]: ');
            if ~isnumeric(energy_list) || isempty(energy_list)
                warning('Invalid input. Using default energy values.');
                energy_list = config.slice_energies;
            end
            
            % Find nearest indices for each energy
            slice_indices = zeros(size(energy_list));
            for i = 1:length(energy_list)
                [~, slice_indices(i)] = min(abs(data.omega_values - energy_list(i)));
                fprintf('Selected energy %.4f eV (slice %d) for QPI\n', energy_list(i), slice_indices(i));
            end
            
            slice_energies = energy_list;
            
        else
            % Index-based selection
            fprintf('Enter QPI slice indices (1-%d):\n', size(data.used, 3));
            
            % Get list of indices
            slice_list = input('Enter list of QPI slice indices [i1, i2, ...]: ');
            if ~isnumeric(slice_list) || isempty(slice_list) || ...
               any(slice_list < 1) || any(slice_list > size(data.used, 3))
                warning('Invalid input. Using default slice indices.');
                slice_list = min(config.slice_indices_direct, size(data.used, 3));
            end
            
            slice_indices = slice_list;
            slice_energies = data.omega_values(slice_list);
        end
    end
end

