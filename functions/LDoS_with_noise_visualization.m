function LDoS_noisy=LDoS_with_noise_visualization(LDoS_result, omega_values, desired_SNR, defect_locations, N)
    % Default single defect location
    if nargin < 4 || isempty(defect_locations)
        defect_locations = [N/2, N/2];
    end

    % Ensure defect_locations is in the correct format
    if size(defect_locations, 1) == 1
        defect_locations = repmat(defect_locations, 1, 1);
    end
    
    gridSize= size(LDoS_result,1);
    % Generate noisy data according to SNR
    LDoS_noisy = addNoise(LDoS_result, desired_SNR);

    % Save the result
    save('LDoS_result_multi_defect_noisy.mat', 'LDoS_noisy', 'omega_values', 'desired_SNR', 'defect_locations');

    % Visualize the 3D grid of LDoS results
    rangeType = 'dynamic'; % 'dynamic' or 'global'
    d3gridDisplay(LDoS_noisy, rangeType);

    % Visualize the LDoS at selected energy levels with defect locations
    gridDisplay(LDoS_noisy, omega_values, defect_locations, gridSize, N);
end
