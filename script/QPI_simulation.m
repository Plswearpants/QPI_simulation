%% Descrption
% This script simulates multi-defect defect QPI signals on a given crystal structure. It is based on the
% TB-model and only consider nearest neighbor hopping term. For more
% details, please consult the supplementary file of the pape:
% (doi:https://doi.org/10.1038/s41467-020-14633-1)
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~Computing_LDoS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Define Parameters
a = 1*10^-9; % lattice constant
t = -0.2; % hopping parameter
E0 = 0; % on-site energy
N = 25; % number of lattice points along one dimension

%% 2. Set up simulating ranges
n = 300; % number of grid points for numerical integration
epsilon = 1e-3; % small imaginary part for numerical stability
gridSize = 51; % number of sampling points along one dimension
omega_values = linspace(-0.5, 0.5, 3); % energy levels

%% 3. Define Defect configuration
simulation_type = 'single'; % Options: 'single' or 'multi'
if strcmp(simulation_type, 'single')
    num_defects = 1;
    defect_energies = -0.02; % single defect energy
    defect_location = [ceil(N/2), ceil(N/2)]; % center position for single defect
else
    num_defects = 3;
    defect_energies = -0.02 * ones(num_defects, 1); % multiple defects with same energy
    defect_location = assignDefectLocations(num_defects, N); % randomly assign defect locations
end

%% 3. Compute LDoS
tic;
disp(['Computing ' simulation_type ' defect LDoS...']);

% Compute LDoS using computeLDoSCore 
[LDoS_result, used_locations] = computeLDoSCore(omega_values, defect_energies, defect_location, ...
    N, a, t, E0, n, epsilon, gridSize);
elapsed_time = toc;
disp(['Computation completed in ' num2str(elapsed_time) ' seconds']);

% Save the result
save_filename = ['LDoS_' simulation_type '_defect.mat'];
save(save_filename, 'LDoS_result', 'used_locations', 'omega_values', 'epsilon', 'n', 'N');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~Loading~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('LDoS_result.mat', 'LDoS_result', 'omega_values', 'epsilon', 'n');
N=50;
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~Visualization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Plot the band structure and constant energy contours

% Define parameters for the dispersion relation
a = 1e-9; % lattice constant
t = -0.2; % hopping parameter
E0 = 0; % on-site energy

% Define k-space grid
k_vals = linspace(-pi/a, pi/a, 100);
[kx, ky] = meshgrid(k_vals, k_vals);

% Compute the energy dispersion
E = E0 - 2 * t * (cos(kx * a) + cos(ky * a));

% Plot the energy dispersion
figure;
surf(kx, ky, E);
colorbar;
title('Energy Dispersion Relation');
xlabel('k_x');
ylabel('k_y');
zlabel('Energy (E)');
grid on;

%% Compute and plot CEC and its convolution at different energies
% Define energy range matching the band structure
omega_vals = linspace(min(E(:)), max(E(:)), 41);

% Initialize arrays to store CEC and convolution results
CEC = zeros(size(kx,1), size(kx,2), length(omega_vals));
CEC_conv = zeros(size(kx,1), size(kx,2), length(omega_vals));

% Compute CEC and convolution for each energy
for i = 1:length(omega_vals)
    % Create constant energy contour with small tolerance
    CEC(:,:,i) = double(abs(E - omega_vals(i)) < 0.05);
    
    % Compute convolution using conv2
    temp_conv = conv2(CEC(:,:,i), CEC(:,:,i), 'same');
    CEC_conv(:,:,i) = temp_conv / max(temp_conv(:)); % Normalize
end

figure;
d3gridDisplay(CEC,'dynamic');
figure;
d3gridDisplay(CEC_conv,'dynamic');

%% Plot DOS at different energy(surface area of the Fermi contour)

% Define parameters for the dispersion relation
a = 1; % lattice constant
t = -0.2; % hopping parameter
E0 = 0; % on-site energy
tic;
% Define k-space grid
k_vals = linspace(-pi, pi, 5000); % Finer grid for better accuracy
[kx, ky] = meshgrid(k_vals, k_vals);

% Calculate the DOS as a function of energy
DOS = zeros(length(omega_values), 1);

for i = 1:length(omega_values)
    E = omega_values(i) + E0;
    % Calculate the dispersion relation
    E_k = E0 - 2 * t * (cos(kx * a) + cos(ky * a));
    % Calculate the DOS proportional to the number of (kx, ky) pairs that satisfy E_k = E
    DOS(i) = sum(abs(E_k - E) < 1e-3, 'all'); % Sum grid points within a small tolerance
end

% Normalize the DOS
DOS = DOS / max(DOS);
toc;
% Plot the DOS as a function of energy
figure;
plot(omega_values, DOS, '-o', 'LineWidth', 2);
title('Density of States (DOS) as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('Normalized DOS');
grid on;

%% Trim the dataset
selected_energy= [1:20,22:41];
LDoS_result=LDoS_result(:,:,selected_energy);
%% Generate QPI from target_LDoS
target_LDoS= LDoS_sim;

QPI_sim= zeros(size(target_LDoS));
for k=1:size(target_LDoS,3)
    QPI_sim(:,:,k)=abs(fftshift(fft2(target_LDoS(:,:,k) - mean(mean(target_LDoS(:,:,k))))));
end
%% Visualize LDoS for a specific energy level(single defect)
% -------------------------------------------
energy_level = 40; % example energy level index

figure()
load('InverseGray','invgray')
map = gray;
[X, Y] = meshgrid(linspace(-N*lattice_spacing/2, N*lattice_spacing/2, 301), ...
                                  linspace(-N*lattice_spacing/2, N*lattice_spacing/2, 301));
imagesc(X(1,:),Y(:,1),LDoS_result(:,:,energy_level));
colormap(map); % Set the colormap to invgray
colorbar;
title(['LDoS at \omega - E0 = ', num2str(omega_values(energy_level)-E0)]);
xlabel('x');
ylabel('y');
hold on;

% Add green dots at lattice positions
% -------------------------------------------
lattice_spacing = a;
[lattice_x, lattice_y] = meshgrid(linspace(-N*lattice_spacing/2, N*lattice_spacing/2, N+1), ...
                                  linspace(-N*lattice_spacing/2, N*lattice_spacing/2, N+1));
%plot(lattice_x(:), lattice_y(:), 'g.', 'MarkerSize', 2);

hold off;

%% Visualize LDoS for a specific energy level(multi-defect) with defect location 
% Visualize the defect locations
N=50;
figure;
plot(used_locations(:,1), used_locations(:,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlim([0, N+1]);
ylim([0, N+1]);
title('Defect Locations');
xlabel('x');
ylabel('y');
grid on;
axis equal;

% Visualize the LDoS at selected energy levels
omega_indices = [1, 2, 35]; % Start, middle, and end
selected_omega_values = omega_values(omega_indices);

figure;
for i = 1:length(omega_indices)
    subplot(1, length(omega_indices), i);
    imagesc(LDoS_result(:,:,omega_indices(i)));
    colorbar;
    title(['LDoS at \omega = ', num2str(selected_omega_values(i))]);
    xlabel('x');
    ylabel('y');
    axis square;
end

sgtitle('LDoS at Selected Energy Levels for Multiple Defects')


%% Plot energy dispersion of the grid

%nos=5;
%figure(3)
figure();
load('InverseGray','invgray')
map = gray;

% Convert LDoS_result to image format using the colormap
for k = 1:size(LDoS_result,3)
    MaxValue = max(LDoS_result(:,:,k), [], 'all');
    MinValue = min(LDoS_result(:,:,k), [], 'all');
    % Using dynamic range for each slice (can be adjusted as needed)
    slicedglob_LDoS(:,:,k,:) = mat2im(LDoS_result(:,:,k), map, [1.2*MinValue 0.8*MaxValue]);
end

imshow3D(slicedglob_LDoS);



%% plot QPI dispersion 
figure();
d3gridDisplay(QPI_sim,'dynamic')

%% Plot 3D QPI dispersion
figure;
[kx, ky] = meshgrid(linspace(-pi, pi, size(QPI_sim,1)), linspace(-pi, pi, size(QPI_sim,2)));

% Plot for a selected energy slice (e.g. middle of the energy range)
energy_slice = ceil(size(QPI_sim,3)/2);
surf(kx, ky, QPI_sim(:,:,energy_slice));
colorbar;
xlabel('k_x');
ylabel('k_y');
zlabel('Intensity');
grid on;


%% Plot combined LDoS & QPI 

omega_values = linspace(-0.5, 0.5, 41);
path = 'Simulated_LDoS&QPI'; % Change your file name here
mkdir(path);
load('InverseGray', 'invgray'); % Ensure 'invgray' colormap is available

% Assign datasets to plot
Grid = LDoS_result; % Ensure 'LDoS_result' variable is loaded
qpi = QPI_sim; % Ensure 'QPI_sim' variable is loaded

for k = 1:length(omega_values)
    j = k;
    f = figure();
    f.Position(3:4) = [2400 1000]; % Increase the figure size for higher resolution
    
    % Make tiles 
    tiledlayout(1, 2);
    
    % Tile 1 
    nexttile;
    imagesc(Grid(:,:,k));
    colormap(invgray); % Use 'invgray' colormap loaded earlier
    title([num2str(1000 * omega_values(k)), ' mV']);
    axis equal tight;
    colorbar;
    
    % Tile 2 
    nexttile;
    imagesc(qpi(:,:,k));
    colormap(gray); % Use 'gray' colormap for QPI
    axis equal tight;
    colorbar;
    
    % Save the image
    fname = sprintf('QPI_sim%03u_(r,E=%0.3fmV).png', j, omega_values(j) * 1000);
    exportgraphics(f, fullfile(path, fname), 'Resolution', 300); % Set the resolution to 300 DPI
    
    % Close the figure to save memory
    close(f);
end

%% Pictures to Video

% Specify the folder where the PNG files are
imageFolder = 'Simulated_LDoS&QPI'; % Ensure the correct path is set

% Get a list of all PNG files in the folder
imageFiles = dir(fullfile(imageFolder, '*.png'));

% Read the first image to get the dimensions
firstImage = imread(fullfile(imageFolder, imageFiles(1).name));
[height, width, ~] = size(firstImage);

% Create a VideoWriter object
outputVideo = VideoWriter(fullfile(imageFolder, 'outputVideo.mp4'), 'MPEG-4');
outputVideo.FrameRate = 8; % Set the frame rate (frames per second)

open(outputVideo);

% Loop through each image file
for ii = 1:length(imageFiles)
    % Construct the full filename
    imgFile = fullfile(imageFolder, imageFiles(ii).name);
    
    % Read the image
    img = imread(imgFile);
    
    % Resize the image to match the dimensions of the first image
    resizedImg = imresize(img, [height width]);
    
    % Write the image to the video
    writeVideo(outputVideo, resizedImg);
end

% Close the VideoWriter object
close(outputVideo);

disp('Video creation complete.');


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~Processing & Analysis~~~~~~~~~~~~~~~~~~~~~~~~~
%% Decay pattern visualization 

% Select an omega index to visualize
omega_index = 13; % Choose the middle omega value for visualization

% Extract the LDoS for the selected omega value
LDoS_selected = LDoS_result_noisy(:,:,omega_index);

% Grid size
[gridSize, ~] = size(LDoS_selected);

% Center index
center_idx = ceil(gridSize / 2);

% Horizontal profile line across the center
horizontal_profile = LDoS_selected(center_idx, :);

% Diagonal profile line across the center
diagonal_profile = diag(LDoS_selected);

% Generate distance array for decay patterns
r_horizontal = abs((1:gridSize) - center_idx);
r_diagonal = sqrt(2) * r_horizontal;

% Calculate decay patterns to avoid singularity at zero and align max values
max_horizontal_profile = max(horizontal_profile);
max_diagonal_profile = max(diagonal_profile);

decay_1_r_horizontal = 1 ./ (r_horizontal + 1) * max_horizontal_profile;
decay_1_r2_horizontal = 1 ./ ((r_horizontal + 1) .^ 2) * max_horizontal_profile;

decay_1_r_diagonal = 1 ./ (r_diagonal + 1) * max_diagonal_profile ;
decay_1_r2_diagonal = 1 ./ ((r_diagonal + 1) .^ 2) * max_diagonal_profile ;

% Normalize the profiles for comparison
horizontal_profile_normalized = horizontal_profile / max_horizontal_profile;
diagonal_profile_normalized = diagonal_profile / max_diagonal_profile;

% Normalize the decay patterns
decay_1_r_horizontal_normalized = decay_1_r_horizontal / max_horizontal_profile;
decay_1_r2_horizontal_normalized = decay_1_r2_horizontal / max_horizontal_profile;

decay_1_r_diagonal_normalized = decay_1_r_diagonal / max_diagonal_profile;
decay_1_r2_diagonal_normalized = decay_1_r2_diagonal / max_diagonal_profile;

% Create a figure for the profiles
figure;

% Plot the horizontal profile
subplot(2, 1, 1);
plot(horizontal_profile_normalized, 'LineWidth', 2);
hold on;
plot(decay_1_r_horizontal_normalized, '--', 'LineWidth', 2);
plot(decay_1_r2_horizontal_normalized, ':', 'LineWidth', 2);
hold off;
title('Horizontal Profile across the Center');
xlabel('Position');
ylabel('Normalized LDoS');
legend('LDoS', '1/(r+1) decay', '1/(r+1)^2 decay');
grid on;

% Plot the diagonal profile
subplot(2, 1, 2);
plot(diagonal_profile_normalized, 'LineWidth', 2);
hold on;
plot(decay_1_r_diagonal_normalized, '--', 'LineWidth', 2);
plot(decay_1_r2_diagonal_normalized, ':', 'LineWidth', 2);
hold off;
title('Diagonal Profile across the Center');
xlabel('Position');
ylabel('Normalized LDoS');
legend('LDoS', '1/(r+1) decay', '1/(r+1)^2 decay');
grid on;

% Add a common title
sgtitle(['Decay Pattern of LDoS for \omega = ', num2str(omega_values(omega_index)), ...
         ', \epsilon = ', num2str(epsilon), ', n = ', num2str(n)]);

%% Zero mean addition and visualization (According to the paper) 
SNR= 6;
selected_energy = omega_values;
% Example usage
LDoS_result_noisy = LDoS_with_noise_visualization(LDoS_result, selected_energy, SNR, [],  N);

%% Zero-mean noise addition (power-based SNR)

% Desired SNR (e.g., 10)
desired_snr = 3;

% Initialize arrays to store signal power and noise power
signal_power = zeros(length(omega_values), 1);
noise_power = zeros(length(omega_values), 1);

% Calculate the signal power for each energy slice and the corresponding noise power
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    signal_power(p) = mean(LDoS_slice(:).^2);
    noise_power(p) = signal_power(p) / desired_snr;
end

% Generate zero-mean noise and add it to the LDoS data
LDoS_result_noisy = zeros(size(LDoS_result));
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    % Generate noise with the calculated noise power
    noise_std = sqrt(noise_power(p));
    noise = noise_std * randn(size(LDoS_slice));
    % Add noise to the LDoS data
    LDoS_result_noisy(:,:,p) = LDoS_slice + noise;
end


% Plot signal power and noise power for verification
figure;
plot(omega_values, signal_power, '-o', 'LineWidth', 2);
hold on;
plot(omega_values, noise_power, '-x', 'LineWidth', 2);
hold off;
title('Signal Power and Noise Power as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('Power');
legend('Signal Power', 'Noise Power');
grid on;
%% Zero-mean noise addition (amplitude-based SNR)

% Desired SNR (e.g., 10)
desired_snr = 0.7; % This is the ratio of RMS signal to RMS noise

% Initialize arrays to store signal RMS and noise RMS
signal_rms = zeros(length(omega_values), 1);
noise_rms = zeros(length(omega_values), 1);

% Calculate the RMS of the LDoS data and the required noise RMS for the desired SNR
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    signal_rms(p) = rms(LDoS_slice(:));
    noise_rms(p) = signal_rms(p) / desired_snr;
end

% Generate zero-mean noise and add it to the LDoS data
LDoS_result_noisy = zeros(size(LDoS_result));
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    % Generate noise with the calculated noise RMS
    noise_std = noise_rms(p);
    noise = noise_std * randn(size(LDoS_slice));
    % Add noise to the LDoS data
    LDoS_result_noisy(:,:,p) = LDoS_slice + noise;
end

% Save the noisy result for future use
save('LDoS_result_noisy.mat', 'LDoS_result_noisy', 'omega_values', 'epsilon', 'n');

% Verify the SNR for each energy slice
snr_values = zeros(length(omega_values), 1);
for p = 1:length(omega_values)
    original_slice = LDoS_result(:,:,p);
    noisy_slice = LDoS_result_noisy(:,:,p);
    noise_slice = noisy_slice - original_slice;
    snr_values(p) = computeSNR(original_slice, noise_slice);
end

% Plot SNR values to verify they are close to the desired SNR
figure;
plot(omega_values, snr_values, '-o', 'LineWidth', 2);
title('SNR of Noisy Data as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('SNR');
grid on;

%% Zero-mean noise addition (frequency-domain SNR)

% Desired SNR (e.g., 10)
desired_snr = 10;

% Initialize arrays to store signal power and noise power
signal_power_spectrum = zeros(length(omega_values), 1);
noise_power_spectrum = zeros(length(omega_values), 1);

% Calculate the Fourier transform of the LDoS data and its power spectrum
LDoS_ft = fftshift(fft2(LDoS_result));

% Calculate the signal power spectrum and the corresponding noise power spectrum
for p = 1:length(omega_values)
    LDoS_ft_slice = LDoS_ft(:,:,p);
    signal_power_spectrum(p) = mean(abs(LDoS_ft_slice(:)).^2);
    noise_power_spectrum(p) = signal_power_spectrum(p) / desired_snr;
end

% Generate noise in the frequency domain and add it to the LDoS data
LDoS_result_noisy = zeros(size(LDoS_result));
for p = 1:length(omega_values)
    LDoS_ft_slice = LDoS_ft(:,:,p);
    % Generate noise in the frequency domain with the calculated noise power
    noise_std = sqrt(noise_power_spectrum(p));
    noise_ft = noise_std * (randn(size(LDoS_ft_slice)) + 1i * randn(size(LDoS_ft_slice)));
    % Add noise in the frequency domain
    LDoS_ft_noisy_slice = LDoS_ft_slice + noise_ft;
    % Transform back to the spatial domain
    LDoS_result_noisy(:,:,p) = real(ifft2(ifftshift(LDoS_ft_noisy_slice)));
end

% Save the noisy result for future use
save('LDoS_result_noisy.mat', 'LDoS_result_noisy', 'omega_values', 'epsilon', 'n');

% Plot signal power spectrum and noise power spectrum for verification
figure;
plot(omega_values, signal_power_spectrum, '-o', 'LineWidth', 2);
hold on;
plot(omega_values, noise_power_spectrum, '-x', 'LineWidth', 2);
hold off;
title('Signal Power Spectrum and Noise Power Spectrum as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('Power Spectrum');
legend('Signal Power Spectrum', 'Noise Power Spectrum');
grid on;
%% Zero-mean noise addition (Peak SNR)

% Desired PSNR (e.g., 30 dB)
desired_psnr = 30;

% Maximum possible pixel value (assuming normalized data in [0, 1])
max_value = 1;

% Initialize an array to store the noise variance for each energy slice
noise_variance = zeros(length(omega_values), 1);

% Calculate the required noise variance for the desired PSNR
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    mse = (max_value^2) / (10^(desired_psnr / 10));
    noise_variance(p) = mse;
end

% Generate zero-mean noise and add it to the LDoS data
LDoS_result_noisy = zeros(size(LDoS_result));
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    % Generate noise with the calculated variance
    noise_std = sqrt(noise_variance(p));
    noise = noise_std * randn(size(LDoS_slice));
    % Add noise to the LDoS data
    LDoS_result_noisy(:,:,p) = LDoS_slice + noise;
end

% Save the noisy result for future use
save('LDoS_result_noisy.mat', 'LDoS_result_noisy', 'omega_values', 'epsilon', 'n');


% Verify the PSNR for each energy slice
psnr_values = zeros(length(omega_values), 1);
for p = 1:length(omega_values)
    original_slice = LDoS_result(:,:,p);
    noisy_slice = LDoS_result_noisy(:,:,p);
    psnr_values(p) = computePSNR(original_slice, noisy_slice, max_value);
end

% Plot the PSNR values to verify they are close to the desired PSNR
figure;
plot(omega_values, psnr_values, '-o', 'LineWidth', 2);
title('PSNR of Noisy Data as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('PSNR (dB)');
grid on;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Tests~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% singal power at different energy levels

% Initialize an array to store the signal power for each energy slice
signal_power = zeros(length(omega_values), 1);

% Compute the signal power for each energy slice
for p = 1:length(omega_values)
    LDoS_slice = LDoS_result(:,:,p);
    signal_power(p) = mean(LDoS_slice(:).^2);
end

% normalize the signal power 
signal_power=signal_power/max(signal_power, [] ,"all");

% compute the signal strength 
signal_strength = sqrt(signal_power);
% Plot the signal power as a function of energy
figure;
plot(omega_values, signal_strength, '-o', 'LineWidth', 2);
title('Signal Power as a Function of Energy');
xlabel('Energy (\omega)');
ylabel('Signal Power');
grid on;

%% T_matrix at different energy levels

% Define parameters
a = 1; % lattice constant
t = -0.2; % hopping parameter
E0 = 0; % on-site energy
Ed = 0; % defect energy
n = 500; % number of points for numerical integration
epsilon = 1e-3; % small imaginary part for numerical stability

% Define omega values
omega_values = linspace(-0.5, 0.5, 41);

% Initialize an array to store T-matrix values
T_values = zeros(size(omega_values));

% Compute the T-matrix values for each omega
for i = 1:length(omega_values)
    omega = omega_values(i);
    T = computeTMatrix(Ed, omega, a, t, E0, n, epsilon);
    T_values(i) = T; % Since T is a scalar, we can directly store it
end

% Plot T-matrix values vs. omega
figure;
plot(omega_values, real(T_values), '-o', 'LineWidth', 2);
hold on;
plot(omega_values, imag(T_values), '-x', 'LineWidth', 2);
hold off;
title('T-matrix Value vs. \omega');
xlabel('Energy (\omega)');
ylabel('T-matrix Value');
legend('Real part', 'Imaginary part');
grid on;

%% Functions dedicated only to this script

% Function to compute RMS
function rms_value = rms(signal)
    rms_value = sqrt(mean(signal.^2));
end

% Function to compute SNR for verification
function snr_value = computeSNR(signal, noise)
    signal_rms = sqrt(mean(signal.^2));
    noise_rms = sqrt(mean(noise.^2));
    snr_value = signal_rms / noise_rms;
end

% Function to compute PSNR for verification
function psnr_value = computePSNR(original, distorted, max_value)
    % Compute the Mean Squared Error (MSE)
    mse = mean((original(:) - distorted(:)).^2);
    
    % Compute PSNR
    if mse == 0
        psnr_value = Inf; % If MSE is zero, PSNR is infinite (no error)
    else
        psnr_value = 10 * log10((max_value^2) / mse);
    end
end