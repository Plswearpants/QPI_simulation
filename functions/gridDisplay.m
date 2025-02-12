function gridDisplay(LDoS_noisy, omega_values, defect_locations, gridSize, N)
    % Select energy levels to visualize
    omega_indices = [1, floor(length(omega_values)/2), length(omega_values)]; % Start, middle, and end
    selected_omega_values = omega_values(omega_indices);

    % Adjust defect locations for plotting
    adjusted_defect_locations = defect_locations * (gridSize / N);

    % Visualize the LDoS at selected energy levels with defect locations
    figure;
    for i = 1:length(omega_indices)
        subplot(1, length(omega_indices), i);
        imagesc(LDoS_noisy(:,:,omega_indices(i)));
        colorbar;
        hold on;
        % Adjust for image coordinates
        scatter(adjusted_defect_locations(:,1), gridSize - adjusted_defect_locations(:,2) + 1, 100, 'r', 'filled');
        hold off;
        title(['Noisy LDoS at \omega = ', num2str(selected_omega_values(i))]);
        xlabel('x');
        ylabel('y');
        axis equal tight;
        set(gca, 'YDir', 'normal'); % Correct the y-axis direction
    end

    sgtitle('Noisy LDoS at Selected Energy Levels for Multiple Defects');
end