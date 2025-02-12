function noise_variance = estimate_noise(image, method)
    % Convert the image to grayscale if it is a color image
    if size(image, 3) == 3
        img_gray = rgb2gray(image);
    else
        img_gray = image;
    end

    % Check the selected method and proceed accordingly
    switch lower(method)
        case 'std'
            % Display the image
            figure;
            imshow(img_gray,[]);
            title('Draw a rectangle to select the region of interest (ROI) of noise. Double-click to confirm.');

            % Let the user draw a rectangle to select the ROI
            h = imrect;
            position = wait(h);  % Wait until the ROI is double-clicked

            % Extract the region of interest (ROI)
            roi = imcrop(img_gray, position);

            % Calculate the variance as the noise level
            noise_variance = var(double(roi(:)));

            % Display the result
            disp(['Estimated noise variance (standard deviation method): ', num2str(noise_variance)]);

            % Close the figure
            close;
            
        case 'wavelet'
            % Convert the image to double precision for wavelet transform
            img_double = double(img_gray);

            % Perform wavelet decomposition using the 'db1' wavelet (Haar wavelet)
            [cA, cH, cV, cD] = dwt2(img_double, 'db1');

            % Calculate the variance of the high-frequency components
            noise_var_H = var(cH(:));
            noise_var_V = var(cV(:));
            noise_var_D = var(cD(:));

            % Combine the variances as a measure of the noise level
            noise_variance = noise_var_H + noise_var_V + noise_var_D;

            % Display the result
            disp(['Estimated noise variance (wavelet method): ', num2str(noise_variance)]);

        otherwise
            error('Invalid method. Choose either "std" or "wavelet".');
    end
end
