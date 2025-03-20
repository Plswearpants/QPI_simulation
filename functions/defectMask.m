function [dIdV_masked] = defectMask(dIdV, midV, V, varargin)
% Description: 
%   Function takes in a 3D dIdV matrix (x,y,V)
%   Also takes voltage vector, and user specified voltage to make mask at
%   Outputs dIdV matrix with Gaussian smooth at all specified defect centres
%   Also outputs FT QPI map of masked dIdV

% Input: 
%   dIdV: dIdV from load3dsall
%   midV: averaged bias array
%   V: the bias voltage at which we pick the defects from
%   varargin: optional input locations [x,y] as an Nx2 matrix
% Output: 
%   dIdV_masked: dIdV after the masking(same as dIdV)
%   QPI_masked: 2D FFT of dIdV_masked

% Parse input arguments
p = inputParser;
addOptional(p, 'locations', [], @(x) isempty(x) || (isnumeric(x) && size(x,2) == 2));
parse(p, varargin{:});
input_locations = p.Results.locations;

%% Set dIdV to average 0

% Remove dIdV offset
for i = 1:size(dIdV,3)
    dIdV(:,:,i) = dIdV(:,:,i) - mean(mean(dIdV(:,:,i)));
end

%%

% Find index closest to user specified voltage
[~, index] = min(abs(V - midV(:)));

% set colour scale limits
clims = ([min(min(dIdV(:,:,index))) 0.1*max(max(dIdV(:,:,index)))]);

% If locations are provided, use them directly; otherwise use interactive mode
if ~isempty(input_locations)
    % Use the provided locations
    x = input_locations(:,1)';
    y = input_locations(:,2)';
    
    % Optionally show the image with marked points for confirmation
    figure;
    imagesc(dIdV(:,:,index), clims);
    colormap("gray");
    pbaspect([1 1 1]);
    hold on;
    scatter(x, y, 50, 'r', 'filled', 'MarkerEdgeColor', 'k');
    title('Defect locations (from input)');
    drawnow;
    pause(1); % Show for a moment
    
    fprintf('Using %d provided defect locations\n', length(x));
else
    % Interactive mode - let user select points
    figure;
    imagesc(dIdV(:,:,index), clims);
    colormap("gray");
    pbaspect([1 1 1]);
    title('Click on defects to mask, then press Enter when done');
    
    i = 1;
    clear x
    clear y
    
    hold on;
    % Until user ends process by typing 'n' apply Gaussian mask at pointer
    % location. Continue by typing 'y'. Creates two vectors with x and y index
    % of gaussian masks
    while 1
        [x(i), y(i)] = ginput(1);
        if isempty(x(i)) % User pressed Enter without clicking
            x = x(1:i-1);
            y = y(1:i-1);
            break;
        end
        h = drawpoint('Position', [x(i) y(i)]);
        resp = input('Do you wish to mask another point? y/n \n', 's');
        i = i + 1;
        fprintf('Selected point %d\n', i-1); % Outputs defect number
        if strcmpi(resp, 'n')
            % user has typed in N or n so break out of the while loop
            break;
        end
    end
    hold off;
end

% Check if any points were selected
if isempty(x)
    warning('No defects selected. Returning original data.');
    dIdV_masked = dIdV;
    return;
end

sigma = 10; % Standard deviation for Gaussian mask
B = zeros(size(dIdV,1), size(dIdV,2), length(x)); % Create masking matrix
dIdV_masked = dIdV; % Create output dIdV (initialize to input data)

fprintf('Creating mask for %d defects\n', length(x));

% Creates a Gaussian mask at each specified point
for k = 1:length(x)
    for i = 1:size(dIdV,1)
        for j = 1:size(dIdV,2)
            B(i,j,k) = 1 - exp(-(j-round(x(k)))^2/(2*sigma^2) - (i-round(y(k)))^2/(2*sigma^2));
        end
    end
end

fprintf('Applying masks to data\n');

% Multiplies the dIdV by the mask, removing defect signatures
for k = 1:length(x)
    dIdV_masked = dIdV_masked .* B(:,:,k);
end

fprintf('Masking complete\n');

end
