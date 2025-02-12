% Load and prepare data
data = LDoS_result;  % 301x301x41 dataset
energy_points = 41;
energy_range = linspace(-20, 20, energy_points);

% Create figure with three subplots (two original + one peak detection visualization)
figure('Position', [100 100 1800 500]);

%% First subplot - Diagonal mask analysis
subplot(1,3,1);

% Create and apply diagonal mask
mask_diag = eye(301);
mask_diag(1:151, 1:151) = 0;  % Zero out left half
masked_data_diag = zeros(size(data));
for t = 1:energy_points
    masked_data_diag(:,:,t) = data(:,:,t) .* mask_diag;
end

% Extract diagonal elements
q_diag = zeros(150, energy_points);
for t = 1:energy_points
    diag_elements = diag(masked_data_diag(:,:,t));
    q_diag(:,t) = diag_elements(152:301);
end

% Calculate envelope and lifetime using findpeaks
lifetime_diag = zeros(1, energy_points);
envelope_diag = zeros(size(q_diag));
threshold = 1/exp(1);

for e = 1:energy_points
    % Find peaks
    [peaks, locs] = findpeaks(q_diag(:,e));
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
imagesc(energy_range, 1:150, q_diag);
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
mask_horiz = zeros(301);
mask_horiz(151,:) = 1;
masked_data_horiz = zeros(size(data));
for t = 1:energy_points
    masked_data_horiz(:,:,t) = data(:,:,t) .* mask_horiz;
end

% Extract horizontal elements
q_horiz = zeros(150, energy_points);
for t = 1:energy_points
    full_line = masked_data_horiz(151,:,t);
    q_horiz(:,t) = full_line(152:301);
end

% Calculate envelope and lifetime using findpeaks
lifetime_horiz = zeros(1, energy_points);
envelope_horiz = zeros(size(q_horiz));

for e = 1:energy_points
    % Find peaks
    [peaks, locs] = findpeaks(q_horiz(:,e));
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
imagesc(energy_range, 1:150, q_horiz);
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

% Choose middle energy slice (e=21)
middle_slice = q_diag(:,11);  % Using diagonal data for visualization
[peaks, locs] = findpeaks(middle_slice);

% Plot the data and peaks
plot(1:150, middle_slice, 'b-', 'LineWidth', 1);
hold on;
plot(locs, peaks, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
% Plot the envelope
if ~isempty(peaks)
    envelope = interp1(locs, peaks, 1:150, 'pchip', 'extrap');
    plot(1:150, envelope, 'r--', 'LineWidth', 2);
end
xlabel('Position');
ylabel('Intensity');
title('Peak Detection (E = 0)');
grid on;
hold off;

