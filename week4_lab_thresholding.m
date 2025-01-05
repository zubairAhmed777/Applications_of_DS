
% Function to calculate pixel intensity distribution
function result = customTabulateAll(A)
    % Flatten the matrix to a single column vector
    A = A(:);

    %disp(A);
    % Find max value 
    maxVal = max(A);
    
    % Preallocate space for the result array
    %result = zeros(maxVal, 3);
    result = zeros(maxVal, 3);

    % Calculate total number of elements
    %totalElements = numel(A);
    
    % Loop over each value to compute count and percentage
    for i = 0:maxVal
        
        % Count occurrences of the value
        count = sum(A == i);
        
        % Calculate percentage
        %percentage = (count / totalElements) * 100;
        
        % Store the value, count, and percentage in the result array
        %result(i, :) = [i, count, percentage];
        result(i+1, :) = [i+1, i, count];
    end
end

% Function to compute within-class variance
function result = withinClassVariance(hist_matrix, threshold, total_values_c)
    weight_sum_bg = 0;
    mean_sum_bg = 0;
    total_values = hist_matrix(end, 1);
    %m_val_count = 0;
    for i=0:threshold-1
        weight_sum_bg = weight_sum_bg + hist_matrix(i+1, 3);
        mean_sum_bg = mean_sum_bg + ( hist_matrix(i+1, 2) * hist_matrix(i+1, 3));
        %m_val_count = m_val_count + 
    end
    %disp(weight_sum_bg);
    weight_bg = weight_sum_bg/total_values_c;
    %disp(weight_bg);
    mean_bg = mean_sum_bg/weight_sum_bg;
    variance_bg = 0;
    for i=0:threshold-1
        temp = hist_matrix(i+1, 2) - mean_bg;
        temp = (temp * temp) * hist_matrix(i+1, 3);
        variance_bg = variance_bg + temp;
    end

    variance_bg = variance_bg/weight_sum_bg;

    weight_sum_fg = 0;
    mean_sum_fg = 0;
    
    %m_val_count = 0;
    for i=threshold:total_values-1
        weight_sum_fg = weight_sum_fg + hist_matrix(i+1, 3);
        mean_sum_fg = mean_sum_fg + ( hist_matrix(i+1, 2) * hist_matrix(i+1, 3));
        %m_val_count = m_val_count + 
    end
    weight_fg = weight_sum_fg/total_values_c;
    mean_fg = mean_sum_fg/weight_sum_fg;
    variance_fg = 0;
    for i=threshold:total_values-1
        temp = hist_matrix(i+1, 2) - mean_fg;
        temp = (temp * temp) * hist_matrix(i+1, 3);
        variance_fg = variance_fg + temp;
    end
    
    % Handle NaN values for empty foreground or background
    if isnan(variance_bg)
        variance_bg = 0;
    end

    if isnan(variance_fg)
        variance_fg = 0;
    end


    variance_fg = variance_fg/weight_sum_fg;
    %{
    disp('weight bg');
    disp(weight_bg);
    disp('var bg: ')
    disp(variance_bg);
    disp('weight fg');
    disp(weight_fg)
    disp('var fg');
    disp(variance_fg);
    %}
    result = (weight_fg * variance_fg) + (weight_bg * variance_bg);
end

% Function to apply thresholding
function result = custThreshold(img, threshold_val)
    result=zeros(size(img));
    [row, col] = size(result);
    for i=1:row
        for j=1:col
            if img(i,j) <= threshold_val
                result(i,j) = 0;
            else
                result(i,j) = 1;
            end
        end
    end
end

% 
%A = [0,0,1,4,4,5; 0,1,3,4,3,4; 1,3,4,2,1,3; 4,4,3,1,0,0; 5,4,2,1,0,0; 5,5,4,3,1,0];

% URL of the raw image file on GitHub
url = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/fp.png';

% Specifying the filename to save the downloaded file
output_filename = 'fp.png';

% Set web options and download the file
opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary

% Download the image file with try-catch block
try
    websave(output_filename, url, opts);
    disp(['File downloaded and saved as ', output_filename]);
catch ME
    disp(['Failed to download the file. Error: ', ME.message]);
end

% Read the image
A = imread(output_filename);

result = customTabulateAll(A);  % Returns a matrix where each row is [value, count, percentage]
all_pixels = (result(:, 2))';
%disp(all_pixels);
variance_all = zeros(size(all_pixels));
%disp(variance_all);
total_values = size(A,1) * size(A, 2);
%disp(total_values);
%Display results
%disp(result(:, :));  % Only showing value and count columns
% Calculate within-class variance for all thresholds
pixel_values_len = size(all_pixels);
for i=1:pixel_values_len(2)
    temp = withinClassVariance(result, all_pixels(i), total_values);
    variance_all(i) = temp;
end
% Find optimal threshold that minimizes variance
[minValue, index] = min(variance_all);
%disp(minValue);
%disp(index);
threshold_val = all_pixels(index);
%disp(threshold_val);
% Apply threshold to the image
img_thresholded = custThreshold(A, threshold_val);
% Displaying the results
subplot(1,2,1);
imshow(output_filename);
title('Original Image');
subplot(1,2,2);
imshow(img_thresholded);
title('Processed Image');