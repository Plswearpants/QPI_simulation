% Load and prepare data
data = -LDoS_result_3; 
energy_points = 3;
energy_range = linspace(-0.5, 0.5, energy_points);

% Create figure with three subplots (two original + one peak detection visualization)
figure('Position', [100 100 1800 500]);

%% First subplot - Diagonal mask analysis
subplot(1,3,1);

% Create and apply diagonal mask
spatial_dim = size(data,1);
center_point = floor((spatial_dim+1)/2);  % Ensure integer with floor
mask_diag = eye(spatial_dim);
mask_diag(1:center_point, 1:center_point) = 0;  % Zero out left half
masked_data_diag = zeros(size(data));
for t = 1:energy_points
    masked_data_diag(:,:,t) = data(:,:,t) .* mask_diag;
end

% Extract diagonal elements
diagonal_length = spatial_dim - center_point;
q_diag = zeros(diagonal_length, energy_points);
for t = 1:energy_points
    diag_elements = diag(masked_data_diag(:,:,t));
    q_diag(:,t) = diag_elements(center_point+1:spatial_dim);
end

% Calculate envelope and lifetime using findpeaks
lifetime_diag = zeros(1, energy_points);
envelope_diag = zeros(size(q_diag));
threshold = 1/exp(1);

for e = 1:energy_points
    % Find peaks, but always include the first point as a peak
    [peaks, locs] = findpeaks(q_diag(:,e));
    
    % Add the first point as a peak if it's not already included
    if isempty(locs) || locs(1) > 1
        peaks = [q_diag(1,e); peaks];
        locs = [1; locs];
    end
    
    % Create envelope using all peaks
    if ~isempty(peaks)
        envelope_diag(:,e) = interp1(locs, peaks, 1:size(q_diag,1), 'pchip', 'extrap');
    end
    
    % Calculate lifetime
    initial_intensity = envelope_diag(1,e);
    target_value = initial_intensity * threshold;
    positions = find(envelope_diag(:,e) < target_value, 1);
    
    if ~isempty(positions) && initial_intensity > 0
        lifetime_diag(e) = positions;
    else
        lifetime_diag(e) = NaN;
    end
end

% Plot diagonal analysis
imagesc(energy_range, 1:diagonal_length, q_diag);
colorbar;
xlabel('Energy');
ylabel('Position');
title('Diagonal Mask: Raw Data');
colormap('jet');
hold on;
plot(energy_range, lifetime_diag, 'w-', 'LineWidth', 2);
hold off;

%% Second subplot - Horizontal mask analysis
subplot(1,3,2);

% Create and apply horizontal mask
mask_horiz = zeros(spatial_dim);
mask_horiz(center_point,:) = 1;
masked_data_horiz = zeros(size(data));
for t = 1:energy_points
    masked_data_horiz(:,:,t) = data(:,:,t) .* mask_horiz;
end

% Extract horizontal elements
horizontal_length = spatial_dim - center_point;
q_horiz = zeros(horizontal_length, energy_points);
for t = 1:energy_points
    full_line = masked_data_horiz(center_point,:,t);
    q_horiz(:,t) = full_line(center_point+1:spatial_dim);
end

% Calculate envelope and lifetime using findpeaks
lifetime_horiz = zeros(1, energy_points);
envelope_horiz = zeros(size(q_horiz));

for e = 1:energy_points
    % Find peaks, but always include the first point as a peak
    [peaks, locs] = findpeaks(q_horiz(:,e));
    
    % Add the first point as a peak if it's not already included
    if isempty(locs) || locs(1) > 1
        peaks = [q_horiz(1,e); peaks];
        locs = [1; locs];
    end
    
    % Create envelope using all peaks
    if ~isempty(peaks)
        envelope_horiz(:,e) = interp1(locs, peaks, 1:size(q_horiz,1), 'pchip', 'extrap');
    end
    
    % Calculate lifetime
    initial_intensity = envelope_horiz(1,e);
    target_value = initial_intensity * threshold;
    positions = find(envelope_horiz(:,e) < target_value, 1);
    
    if ~isempty(positions) && initial_intensity > 0
        lifetime_horiz(e) = positions;
    else
        lifetime_horiz(e) = NaN;
    end
end

% Plot horizontal analysis
imagesc(energy_range, 1:horizontal_length, q_horiz);
colorbar;
xlabel('Energy');
ylabel('Position');
title('Horizontal Mask: Raw Data');
colormap('jet');
hold on;
plot(energy_range, lifetime_horiz, 'w-', 'LineWidth', 2);
hold off;

%% Third subplot - Peak Detection Visualization
subplot(1,3,3);

% Choose middle energy slice
middle_energy_index = ceil(energy_points/2);
middle_slice = q_diag(:,middle_energy_index);  % Using diagonal data for visualization

% Find peaks, but always include the first point as a peak
[peaks, locs] = findpeaks(middle_slice);

% Add the first point as a peak if it's not already included
if isempty(locs) || locs(1) > 1
    peaks = [middle_slice(1); peaks];
    locs = [1; locs];
end

% Plot the data and peaks
plot(1:diagonal_length, middle_slice, 'b-', 'LineWidth', 1);
hold on;
plot(locs, peaks, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Plot the envelope
if ~isempty(peaks)
    envelope = interp1(locs, peaks, 1:diagonal_length, 'pchip', 'extrap');
    plot(1:diagonal_length, envelope, 'r--', 'LineWidth', 2);
    
    % Mark the 1/e decay point if it exists
    initial_intensity = envelope(1);
    target_value = initial_intensity * threshold;
    decay_point = find(envelope < target_value, 1);
    
    if ~isempty(decay_point)
        plot(decay_point, envelope(decay_point), 'mo', 'MarkerSize', 10, 'LineWidth', 2);
        plot([decay_point, decay_point], [0, initial_intensity], 'm:', 'LineWidth', 1.5);
        text(decay_point + 5, envelope(decay_point), ['Lifetime = ', num2str(decay_point)], ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end

xlabel('Position');
ylabel('Intensity');
title(['Peak Detection (E = ', num2str(energy_range(middle_energy_index)), ')']);
grid on;
hold off;

%% New section: Plot diagonal and horizontal cuts for all energy slices
figure('Position', [100 100 1800 800]);

% Create subplot grid based on number of energy points
num_rows = 2;  % One row for diagonal, one for horizontal
num_cols = energy_points;

% Plot diagonal cuts
for e = 1:energy_points
    subplot(num_rows, num_cols, e);
    plot(1:diagonal_length, q_diag(:,e), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Add envelope if it exists
    if any(envelope_diag(:,e) > 0)
        plot(1:diagonal_length, envelope_diag(:,e), 'r--', 'LineWidth', 1.5);
        
        % Mark the 1/e decay point if it exists
        if ~isnan(lifetime_diag(e))
            decay_point = lifetime_diag(e);
            decay_value = envelope_diag(decay_point, e);
            plot(decay_point, decay_value, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            
            % Add vertical line at decay point
            plot([decay_point, decay_point], [0, envelope_diag(1,e)], 'r:', 'LineWidth', 1);
        end
    end
    
    title(['Diagonal Cut: E = ', num2str(energy_range(e))]);
    xlabel('Position from Center');
    ylabel('Intensity');
    grid on;
    
    % Add lifetime value as text if available
    if ~isnan(lifetime_diag(e))
        text(0.05, 0.95, ['Lifetime = ', num2str(lifetime_diag(e))], ...
            'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
            'BackgroundColor', [1 1 1 0.7]);
    end
    
    hold off;
end

% Plot horizontal cuts
for e = 1:energy_points
    subplot(num_rows, num_cols, e + num_cols);
    plot(1:horizontal_length, q_horiz(:,e), 'b-', 'LineWidth', 1.5);
    hold on;
    
    % Add envelope if it exists
    if any(envelope_horiz(:,e) > 0)
        plot(1:horizontal_length, envelope_horiz(:,e), 'r--', 'LineWidth', 1.5);
        
        % Mark the 1/e decay point if it exists
        if ~isnan(lifetime_horiz(e))
            decay_point = lifetime_horiz(e);
            decay_value = envelope_horiz(decay_point, e);
            plot(decay_point, decay_value, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            
            % Add vertical line at decay point
            plot([decay_point, decay_point], [0, envelope_horiz(1,e)], 'r:', 'LineWidth', 1);
        end
    end
    
    title(['Horizontal Cut: E = ', num2str(energy_range(e))]);
    xlabel('Position from Center');
    ylabel('Intensity');
    grid on;
    
    % Add lifetime value as text if available
    if ~isnan(lifetime_horiz(e))
        text(0.05, 0.95, ['Lifetime = ', num2str(lifetime_horiz(e))], ...
            'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold', ...
            'BackgroundColor', [1 1 1 0.7]);
    end
    
    hold off;
end

% Adjust spacing between subplots
sgtitle('Diagonal and Horizontal Cuts at All Energy Values', 'FontSize', 14);
