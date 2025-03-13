function compare_QP_lifetime_datasets(datasets, dataset_names, energy_range)
% COMPARE_QP_LIFETIME_DATASETS Compares multiple datasets of LDOS results
%
% Inputs:
%   datasets - Cell array of LDOS datasets to compare
%   dataset_names - Cell array of names for each dataset (for legend)
%   energy_range - Vector of energy values
%
% Example:
%   compare_QP_lifetime_datasets({LDoS_result1, LDoS_result2}, {'Dataset 1', 'Dataset 2'}, linspace(-0.5, 0.5, 3))

% Check inputs
if ~iscell(datasets) || ~iscell(dataset_names)
    error('Datasets and dataset_names must be cell arrays');
end

num_datasets = length(datasets);
if length(dataset_names) ~= num_datasets
    error('Number of dataset names must match number of datasets');
end

energy_points = length(energy_range);

% Define colors for each dataset
colors = lines(num_datasets);

% Create figure for comparison plots
figure('Position', [100 100 1800 800]);

% Create subplot grid
num_rows = 2;  % One row for diagonal, one for horizontal
num_cols = energy_points;

% Initialize arrays to store lifetime values for all datasets
lifetime_diag_all = zeros(num_datasets, energy_points);
lifetime_horiz_all = zeros(num_datasets, energy_points);

% Process each dataset
for d = 1:num_datasets
    % Get current dataset
    data = -datasets{d};
    
    % Create and apply masks
    spatial_dim = size(data, 1);
    center_point = floor((spatial_dim+1)/2);
    
    % Diagonal mask
    mask_diag = eye(spatial_dim);
    mask_diag(1:center_point, 1:center_point) = 0;
    masked_data_diag = zeros(size(data));
    for t = 1:energy_points
        masked_data_diag(:,:,t) = data(:,:,t) .* mask_diag;
    end
    
    % Horizontal mask
    mask_horiz = zeros(spatial_dim);
    mask_horiz(center_point,:) = 1;
    masked_data_horiz = zeros(size(data));
    for t = 1:energy_points
        masked_data_horiz(:,:,t) = data(:,:,t) .* mask_horiz;
    end
    
    % Extract diagonal and horizontal elements
    diagonal_length = spatial_dim - center_point;
    horizontal_length = spatial_dim - center_point;
    
    q_diag = zeros(diagonal_length, energy_points);
    q_horiz = zeros(horizontal_length, energy_points);
    
    for t = 1:energy_points
        diag_elements = diag(masked_data_diag(:,:,t));
        q_diag(:,t) = diag_elements(center_point+1:spatial_dim);
        
        full_line = masked_data_horiz(center_point,:,t);
        q_horiz(:,t) = full_line(center_point+1:spatial_dim);
    end
    
    % Calculate envelopes and lifetimes
    envelope_diag = zeros(size(q_diag));
    envelope_horiz = zeros(size(q_horiz));
    threshold = 1/exp(1);
    
    for e = 1:energy_points
        % Diagonal direction
        [peaks_diag, locs_diag] = findpeaks(q_diag(:,e));
        if isempty(locs_diag) || locs_diag(1) > 1
            peaks_diag = [q_diag(1,e); peaks_diag];
            locs_diag = [1; locs_diag];
        end
        
        if ~isempty(peaks_diag)
            envelope_diag(:,e) = interp1(locs_diag, peaks_diag, 1:size(q_diag,1), 'pchip', 'extrap');
            
            initial_intensity = envelope_diag(1,e);
            target_value = initial_intensity * threshold;
            positions = find(envelope_diag(:,e) < target_value, 1);
            
            if ~isempty(positions) && initial_intensity > 0
                lifetime_diag_all(d,e) = positions;
            else
                lifetime_diag_all(d,e) = NaN;
            end
        end
        
        % Horizontal direction
        [peaks_horiz, locs_horiz] = findpeaks(q_horiz(:,e));
        if isempty(locs_horiz) || locs_horiz(1) > 1
            peaks_horiz = [q_horiz(1,e); peaks_horiz];
            locs_horiz = [1; locs_horiz];
        end
        
        if ~isempty(peaks_horiz)
            envelope_horiz(:,e) = interp1(locs_horiz, peaks_horiz, 1:size(q_horiz,1), 'pchip', 'extrap');
            
            initial_intensity = envelope_horiz(1,e);
            target_value = initial_intensity * threshold;
            positions = find(envelope_horiz(:,e) < target_value, 1);
            
            if ~isempty(positions) && initial_intensity > 0
                lifetime_horiz_all(d,e) = positions;
            else
                lifetime_horiz_all(d,e) = NaN;
            end
        end
        
        % Plot diagonal cuts for this energy
        subplot(num_rows, num_cols, e);
        hold on;
        
        % Plot raw data
        plot(1:diagonal_length, q_diag(:,e), 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '-');
        
        % Plot envelope
        if any(envelope_diag(:,e) > 0)
            plot(1:diagonal_length, envelope_diag(:,e), 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '--');
            
            % Mark the 1/e decay point if it exists
            if ~isnan(lifetime_diag_all(d,e))
                decay_point = lifetime_diag_all(d,e);
                decay_value = envelope_diag(decay_point, e);
                plot(decay_point, decay_value, 'o', 'Color', colors(d,:), 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
        
        % Plot horizontal cuts for this energy
        subplot(num_rows, num_cols, e + num_cols);
        hold on;
        
        % Plot raw data
        plot(1:horizontal_length, q_horiz(:,e), 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '-');
        
        % Plot envelope
        if any(envelope_horiz(:,e) > 0)
            plot(1:horizontal_length, envelope_horiz(:,e), 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '--');
            
            % Mark the 1/e decay point if it exists
            if ~isnan(lifetime_horiz_all(d,e))
                decay_point = lifetime_horiz_all(d,e);
                decay_value = envelope_horiz(decay_point, e);
                plot(decay_point, decay_value, 'o', 'Color', colors(d,:), 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
    end
end

% Finalize plots with titles, labels, and legends
for e = 1:energy_points
    % Diagonal subplot
    subplot(num_rows, num_cols, e);
    title(['Diagonal Cut: E = ', num2str(energy_range(e))]);
    xlabel('Position from Center');
    ylabel('Normalized Intensity');
    grid on;
    
    % Create normalized plots for diagonal cuts
    for d = 1:num_datasets
        % Get the data for this dataset and energy
        current_data = q_diag(:,e);
        current_envelope = envelope_diag(:,e);
        
        % Normalize by the first point if it's not zero
        if current_data(1) ~= 0
            normalized_data = current_data / current_data(1);
            % Clear the current plot and replot with normalized data
            plot(1:diagonal_length, normalized_data, 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '-');
            hold on;
        end
        
        % Normalize and replot the envelope if it exists
        if any(current_envelope > 0) && current_envelope(1) > 0
            normalized_envelope = current_envelope / current_envelope(1);
            plot(1:diagonal_length, normalized_envelope, 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '--');
            
            % Mark the 1/e decay point if it exists
            if ~isnan(lifetime_diag_all(d,e))
                decay_point = lifetime_diag_all(d,e);
                decay_value = normalized_envelope(decay_point);
                plot(decay_point, decay_value, 'o', 'Color', colors(d,:), 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
    end
    
    % Add legend with proper handling of NaN values
    legend_entries = cell(1, num_datasets);
    for d = 1:num_datasets
        if isnan(lifetime_diag_all(d,e))
            lifetime_str = 'N/A';
        else
            lifetime_str = num2str(lifetime_diag_all(d,e), '%.1f');
        end
        legend_entries{d} = [dataset_names{d}, ' (τ=', lifetime_str, ')'];
    end
    %legend(legend_entries, 'Location', 'best');
    
    % Add reference line at 1/e
    yline(1/exp(1), 'k--', '1/e', 'LineWidth', 1, 'Alpha', 0.5);
    ylim([0, 1.1]);
    
    % Horizontal subplot
    subplot(num_rows, num_cols, e + num_cols);
    title(['Horizontal Cut: E = ', num2str(energy_range(e))]);
    xlabel('Position from Center');
    ylabel('Normalized Intensity');
    grid on;
    
    % Create normalized plots for horizontal cuts
    for d = 1:num_datasets
        % Get the data for this dataset and energy
        current_data = q_horiz(:,e);
        current_envelope = envelope_horiz(:,e);
        
        % Normalize by the first point if it's not zero
        if current_data(1) ~= 0
            normalized_data = current_data / current_data(1);
            % Clear the current plot and replot with normalized data
            plot(1:horizontal_length, normalized_data, 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '-');
            hold on;
        end
        
        % Normalize and replot the envelope if it exists
        if any(current_envelope > 0) && current_envelope(1) > 0
            normalized_envelope = current_envelope / current_envelope(1);
            plot(1:horizontal_length, normalized_envelope, 'Color', colors(d,:), 'LineWidth', 1.5, 'LineStyle', '--');
            
            % Mark the 1/e decay point if it exists
            if ~isnan(lifetime_horiz_all(d,e))
                decay_point = lifetime_horiz_all(d,e);
                decay_value = normalized_envelope(decay_point);
                plot(decay_point, decay_value, 'o', 'Color', colors(d,:), 'MarkerSize', 8, 'LineWidth', 2);
            end
        end
    end
    
    % Add legend with proper handling of NaN values
    legend_entries = cell(1, num_datasets);
    for d = 1:num_datasets
        if isnan(lifetime_horiz_all(d,e))
            lifetime_str = 'N/A';
        else
            lifetime_str = num2str(lifetime_horiz_all(d,e), '%.1f');
        end
        legend_entries{d} = [dataset_names{d}, ' (τ=', lifetime_str, ')'];
    end
    legend(legend_entries, 'Location', 'best');
    
    % Add reference line at 1/e
    yline(1/exp(1), 'k--', '1/e', 'LineWidth', 1, 'Alpha', 0.5);
    ylim([0, 1.1]);
end

% Add overall title
sgtitle('Comparison of Diagonal and Horizontal Cuts Across Datasets', 'FontSize', 14);

% Create a summary figure for lifetime comparison
figure('Position', [100 100 900 400]);

% Plot diagonal lifetimes
subplot(1, 2, 1);
for d = 1:num_datasets
    plot(energy_range, lifetime_diag_all(d,:), 'o-', 'Color', colors(d,:), 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
end
title('Diagonal Lifetime vs Energy');
xlabel('Energy');
ylabel('Lifetime (position units)');
grid on;
%legend(dataset_names, 'Location', 'best');

% Plot horizontal lifetimes
subplot(1, 2, 2);
for d = 1:num_datasets
    plot(energy_range, lifetime_horiz_all(d,:), 'o-', 'Color', colors(d,:), 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
end
title('Horizontal Lifetime vs Energy');
xlabel('Energy');
ylabel('Lifetime (position units)');
grid on;
%legend(dataset_names, 'Location', 'best');

end
