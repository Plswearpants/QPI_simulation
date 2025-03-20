function [result, freq_diff] = fft_subtract(data1, data2, varargin)
% FFT_SUBTRACT - Subtract two matrices in the frequency domain
%
% Description:
%   Takes two n×n matrices, computes their FFTs (preserving both real and 
%   imaginary parts), subtracts them in the frequency domain, and performs 
%   an inverse FFT to return to the spatial domain. Also visualizes the FFT maps.
%
% Usage:
%   [result, freq_diff] = fft_subtract(data1, data2)
%   [result, freq_diff] = fft_subtract(data1, data2, options)
%
% Inputs:
%   data1, data2  - Input matrices of the same size (n×n)
%   options       - Name-value pairs:
%     'center'      - Whether to center the FFT (default: true)
%     'normalize'   - How to normalize: 'none', 'mean', 'minmax', 'zscore', 'energy' (default: 'energy')
%     'window'      - Apply windowing to reduce edge effects: 'none', 'hann', 'hamming', 'tukey' (default: 'none')
%     'plot'        - Whether to plot the FFT maps (default: true)
%     'log_scale'   - Whether to use log scale for FFT visualization (default: true)
%
% Outputs:
%   result    - The result of inverse FFT after subtraction (spatial domain)
%   freq_diff - The difference in the frequency domain before inverse FFT
%
% Example:
%   [spatial_result, freq_diff] = fft_subtract(matrix1, matrix2, 'normalize', 'energy');

% Parse input arguments
p = inputParser;
addRequired(p, 'data1', @isnumeric);
addRequired(p, 'data2', @isnumeric);
addParameter(p, 'center', true, @islogical);
addParameter(p, 'normalize', 'energy', @(x) ismember(x, {'none', 'mean', 'minmax', 'zscore', 'energy'}));
addParameter(p, 'window', 'none', @(x) ismember(x, {'none', 'hann', 'hamming', 'tukey'}));
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'log_scale', true, @islogical);
parse(p, data1, data2, varargin{:});

center = p.Results.center;
normalize_method = p.Results.normalize;
window_type = p.Results.window;
do_plot = p.Results.plot;
log_scale = p.Results.log_scale;

% Check if matrices are the same size
if ~isequal(size(data1), size(data2))
    error('Input matrices must have the same dimensions');
end

% Check if matrices are square
[n, m] = size(data1);
if n ~= m
    error('Input matrices must be square (n×n)');
end

% Store original data for plotting
orig_data1 = data1;
orig_data2 = data2;

% Normalize data based on selected method
switch normalize_method
    case 'none'
        % No normalization
        fprintf('No normalization applied.\n');
    case 'mean'
        % Subtract mean (zero mean)
        data1 = data1 - mean(data1(:));
        data2 = data2 - mean(data2(:));
        fprintf('Normalization: Subtracted mean from both datasets.\n');
    case 'minmax'
        % Scale to [0,1] range
        data1 = (data1 - min(data1(:))) / (max(data1(:)) - min(data1(:)));
        data2 = (data2 - min(data2(:))) / (max(data2(:)) - min(data2(:)));
        fprintf('Normalization: Scaled both datasets to [0,1] range.\n');
    case 'zscore'
        % Z-score normalization (zero mean, unit variance)
        data1 = (data1 - mean(data1(:))) / std(data1(:));
        data2 = (data2 - mean(data2(:))) / std(data2(:));
        fprintf('Normalization: Applied Z-score normalization to both datasets.\n');
    case 'energy'
        % Normalize by total energy (sum of squares)
        data1 = data1 / sqrt(sum(abs(data1(:)).^2));
        data2 = data2 / sqrt(sum(abs(data2(:)).^2));
        fprintf('Normalization: Scaled both datasets to have unit energy.\n');
end

% Apply window function if requested to reduce edge effects
if ~strcmp(window_type, 'none')
    window = create_2d_window(n, window_type);
    data1 = data1 .* window;
    data2 = data2 .* window;
    fprintf('Applied %s window to reduce edge effects.\n', window_type);
end

% Compute FFT of both matrices
fft_data1 = fft2(data1);
fft_data2 = fft2(data2);

% Apply FFT shifting if centering is requested
if center
    fft_data1_centered = fftshift(fft_data1);
    fft_data2_centered = fftshift(fft_data2);
else
    fft_data1_centered = fft_data1;
    fft_data2_centered = fft_data2;
end

% Plot the original data and FFT maps if requested
if do_plot
    % Create figure with 2 rows and 3 columns
    figure('Position', [100, 100, 1000, 800]);
    
    % Plot original data
    subplot(2, 3, 1);
    imagesc(orig_data1);
    axis square;
    colormap('gray');
    colorbar;
    title('Original Data 1');
    xlabel('X');
    ylabel('Y');
    
    subplot(2, 3, 2);
    imagesc(orig_data2);
    axis square;
    colormap('gray');
    colorbar;
    title('Original Data 2');
    xlabel('X');
    ylabel('Y');
    
    subplot(2, 3, 3);
    imagesc(orig_data1 - orig_data2);
    axis square;
    colormap('gray');
    colorbar;
    title('Original Difference');
    xlabel('X');
    ylabel('Y');
    
    % Prepare FFT data for visualization
    if log_scale
        % Use log scale to enhance visibility of features
        fft1_vis = log(1 + abs(fft_data1_centered));
        fft2_vis = log(1 + abs(fft_data2_centered));
        colormap_title = 'Log(1+|FFT|)';
    else
        % Use linear scale
        fft1_vis = abs(fft_data1_centered);
        fft2_vis = abs(fft_data2_centered);
        colormap_title = '|FFT|';
    end
    
    % Determine common colormap scale for consistent visualization
    cmin = min(min(fft1_vis(:)), min(fft2_vis(:)));
    cmax = max(max(fft1_vis(:)), max(fft2_vis(:)));
    
    % Plot first FFT
    subplot(2, 3, 4);
    imagesc(fft1_vis);
    axis square;
    colormap('jet');
    colorbar;
    caxis([cmin, cmax]);
    title(['FFT of Data 1 (', colormap_title, ')']);
    if center
        xlabel('Frequency (centered)');
        ylabel('Frequency (centered)');
    else
        xlabel('Frequency');
        ylabel('Frequency');
    end
    
    % Plot second FFT
    subplot(2, 3, 5);
    imagesc(fft2_vis);
    axis square;
    colormap('jet');
    colorbar;
    caxis([cmin, cmax]);
    title(['FFT of Data 2 (', colormap_title, ')']);
    if center
        xlabel('Frequency (centered)');
        ylabel('Frequency (centered)');
    else
        xlabel('Frequency');
        ylabel('Frequency');
    end
    
    % Calculate and plot difference
    fft_diff_vis = abs(fft_data1_centered - fft_data2_centered);
    if log_scale
        fft_diff_vis = log(1 + fft_diff_vis);
    end
    
    subplot(2, 3, 6);
    imagesc(fft_diff_vis);
    axis square;
    colormap('jet');
    colorbar;
    title(['FFT Difference (', colormap_title, ')']);
    if center
        xlabel('Frequency (centered)');
        ylabel('Frequency (centered)');
    else
        xlabel('Frequency');
        ylabel('Frequency');
    end
    
    % Add title with normalization information
    norm_title = sprintf('FFT Subtraction with %s Normalization', normalize_method);
    if ~strcmp(window_type, 'none')
        norm_title = [norm_title, sprintf(' and %s Window', window_type)];
    end
    sgtitle(norm_title);
end

% Subtract in frequency domain
freq_diff = fft_data1_centered - fft_data2_centered;

% Perform inverse FFT to get back to spatial domain
if center
    freq_diff_unshifted = ifftshift(freq_diff);
    result = ifft2(freq_diff_unshifted);
else
    result = ifft2(freq_diff);
end

% The result might have small imaginary components due to numerical precision
% If they're very small, we can get rid of them
if max(abs(imag(result(:)))) < 1e-10 * max(abs(real(result(:))))
    result = real(result);
end

end

% Helper function to create window functions
function window = create_2d_window(n, type)
    % Create a 2D window function of specified type
    switch type
        case 'hann'
            win1d = hann(n);
        case 'hamming'
            win1d = hamming(n);
        case 'tukey'
            win1d = tukeywin(n, 0.5); % 0.5 is the taper parameter
        otherwise
            win1d = ones(n, 1); % default to rectangular (no window)
    end
    
    % Create 2D window from outer product
    window = win1d * win1d';
end 