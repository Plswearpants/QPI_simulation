%% Step 1: Define Parameters
a = 1*10^-9; % lattice constant
t = -0.2; % hopping parameter
E0 = 0; % on-site energy
Ed = -5; % defect energy
N = 7; % number of lattice points along one dimension
gridSize = 31; % size of the grid for LDoS calculation

% Grid and Defect Position
[X1, X2] = meshgrid(linspace(-N*a/2, N*a/2, gridSize), linspace(-N*a/2, N*a/2, gridSize));
X = cat(3, X1, X2); % Location vector on the grid
Xd = [0, 0]; % Defect position vector

%% Test script for LDoS convergence with different n values

% Values of n to test for convergence
n_values = [500, 1000, 1500];
omega_values = linspace(-0.5, 0.5, 3); % omega values to test

% Initialize a cell array to store results
LDoS_results = cell(length(n_values), length(omega_values));

% Loop over different omega values
for j = 1:length(omega_values)
    omega = omega_values(j);
    disp(['Testing omega = ', num2str(omega)]);
    
    % Loop over different n values
    parfor idx = 1:length(n_values)
        n = n_values(idx);
        disp(['  Computing LDoS for n = ', num2str(n)]);
        
        % Compute LDoS with the current n value and omega value
        LDoS_results{idx, j} = ComputeLDoS(X, omega, a, t, E0, Ed, Xd, n);
    end
end

%% Testing different epsilon value 

% Values of epsilon to test for convergence
epsilon_values = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
omega_values = linspace(-0.5, 0.5, 11); % omega values to test
n = 500; % Number of points for numerical integration

% Initialize a cell array to store results
LDoS_epsilontest = cell(length(epsilon_values), 1);

% Loop over different omega values
for j = 1:length(omega_values)
    omega = omega_values(j);
    disp(['Testing omega = ', num2str(omega)]);
    
    for idx = 1:length(epsilon_values)
        epsilon = epsilon_values(idx);
        disp(['Computing LDoS for epsilon = ', num2str(epsilon)]);
        
        % Compute LDoS with the current epsilon value
        LDoS_epsilontest{idx,j} = ComputeLDoS(X, omega, a, t, E0, Ed, Xd, n, epsilon);
    end
end

%% Testing different epsilon value & n values

% Values of epsilon & n to test for convergence
omega_values = linspace(-0.5,0.5,4); % omega values to test
epsilon_values = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
n_values = [500,1000,2000,3000]; % Number of points for numerical integration

% Initialize a cell array to store results
LDoS_epsilon_n_test = cell(length(omega_values), length(epsilon_values), length(n_values));

% Loop over different omega values
for p = 1:length(omega_values)
    omega=omega_values(p);
    disp(['Testing omega = ', num2str(omega)]);
    for j = 1:length(epsilon_values)
        epsilon = epsilon_values(j);
        disp(['  Testing epsilon = ', num2str(epsilon)]);
    
        for idx = 1:length(n_values)
            n = n_values(idx);
            disp(['     Computing LDoS for n = ', num2str(n)]);
         
            % Compute LDoS with the current epsilon value
            LDoS_epsilon_n_test{idx,j,p} = ComputeLDoS(X, omega, a, t, E0, Ed, Xd, n, epsilon);
        end
    end
end
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~visualization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Visualization script for LDoS convergence with different epsilon and n values

% Visualization script for LDoS convergence with different epsilon and n values

% Load the colormap (assuming 'invgray' is available)
load('InverseGray', 'invgray');

% Grid and Defect Position (assuming these are still in the workspace)
[X1, X2] = meshgrid(linspace(-N*a/2, N*a/2, gridSize), linspace(-N*a/2, N*a/2, gridSize));

% Ask the user to select an omega value
omega_selected = input('Enter the value of omega to visualize: ');

% Find the index of the selected omega value
[~, p] = min(abs(omega_values - omega_selected));

if isempty(p)
    error('The selected omega value is not within the predefined omega_values.');
end

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 1600, 1000]); % Set the figure size for better visualization

% Loop over different epsilon and n values for the selected omega
for j = 1:length(n_values)
    for idx = 1:length(epsilon_values)
        subplot(length(n_values), length(epsilon_values), (j-1)*length(epsilon_values) + idx);
        LDoS_reshaped = reshape(LDoS_epsilon_n_test{j, idx, p}, gridSize, gridSize); % Reshape to 2D grid
        imagesc(X1(1,:), X2(:,1), LDoS_reshaped);
        %colormap(invgray); % Set the colormap to invgray
        colorbar;
        
        % Titles and labels
        if j == 1
            title(['\epsilon = ', num2str(epsilon_values(idx))]);
        end
        if idx == 1
            ylabel(['n = ', num2str(n_values(j))]);
        end
        %xlabel('x');
        %ylabel('y');
    end
end

% Add a common title
sgtitle(['LDoS Convergence for \omega = ', num2str(omega_selected)]);


%% Plot the results to compare LDoS for different epsilon values
Max=max(cell2mat(LDoS_epsilontest),[],'all');
Min=min(cell2mat(LDoS_epsilontest),[],'all');
figure;
omega_values = linspace(-0.5, 0.5, 11); % omega values to test
omega_num=9;
for idx = 1:length(epsilon_values)
    LDoS_reshaped = reshape(LDoS_epsilontest{idx,omega_num}, gridSize, gridSize); % Example energy level index 1
    subplot(2, 3, idx);
    imagesc(X1(1,:), X2(:,1), LDoS_reshaped);
    %clim([Min,Max])
    colorbar;
    title(['\epsilon = ', num2str(epsilon_values(idx)), '\omega =', num2str(omega_values(omega_num))]);
    xlabel('x');
    ylabel('y');
end
%% Plot the results to compare LDoS for different n values and omega values
figure;
for j = 1:length(omega_values)
    for idx = 1:length(n_values)
        subplot(length(omega_values), length(n_values), (j-1)*length(n_values) + idx);
        LDoS_reshaped = reshape(LDoS_results{idx, j}, gridSize, gridSize); % Reshape to 2D grid
        imagesc(X1(1,:), X2(:,1), LDoS_reshaped);
        colorbar;
        title(['n = ', num2str(n_values(idx)), ', \omega = ', num2str(omega_values(j))]);
        xlabel('x');
        ylabel('y');
    end
end
%% Plot the results to compare LDoS for different n values(global colorscale)
figure;
omega_num=33;
omega_values = linspace(-0.5, 0.5, 41);
n_values = [500, 1000, 1500, 2000, 2500, 3000];

Max=max(cell2mat(LDoS_results(:,omega_num)),[],'all');
Min=min(cell2mat(LDoS_results(:,omega_num)),[],'all');
figure();
for idx = 1:length(n_values)
    LDoS_reshaped = reshape(LDoS_results{idx,omega_num}(:,:,1), gridSize, gridSize); % energy level index omega_num 
    subplot(2, 3, idx);
    imagesc(X1(1,:), X2(:,1), LDoS_reshaped);
    %clim([Min,Max])
    colorbar;
    title(['n = ', num2str(n_values(idx)), ', \omega = ', num2str(omega_values(omega_num))]);
    xlabel('x');
    ylabel('y');
end

%% Plot energy dispersion of the grid
%nos=5;
%figure(3)
figure();
load('InverseGray','invgray')
map = invgray;

iteration_num=6;

for k=1:gridSize
    MaxValue = max(LDoS_results{iteration_num,k},[],'all');
    MinValue = min(LDoS_results{iteration_num,k},[],'all');
    slicedglob_LDoS(:,:,k,:) = mat2im(LDoS_results{iteration_num,k}, map, [1.2*MinValue 0.8*MaxValue]);
end
imshow3D(slicedglob_LDoS)
