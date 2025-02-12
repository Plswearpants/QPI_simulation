function LDoS_noisy = addNoise(LDoS_result, desired_SNR)
    % Compute maximum LDoS value for each layer
    max_LDoS_layer = max(max(LDoS_result, [], 1), [], 2); % Result is 1x1xK
    max_LDoS_layer = reshape(max_LDoS_layer, [1, 1, size(LDoS_result, 3)]); % Reshape to 1x1xK

    % Compute the noise variance based on the desired SNR for each layer
    eta = (max_LDoS_layer / desired_SNR).^2;

    % Generate zero-mean Gaussian noise with variance eta for each layer
    noise = sqrt(eta) .* randn(size(LDoS_result));

    % Add noise to the LDoS results
    LDoS_noisy = LDoS_result + noise;
end