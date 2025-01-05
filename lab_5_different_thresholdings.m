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

% Entropy based thresholding
function pt = calculate_pt(matrix)
    % Flatten the matrix into a single vector
    flattened = matrix(:);
    total_pixels = numel(flattened); % Total number of pixels
    result = customTabulateAll(matrix);
    vals = result(:,3)';
    pt = vals/total_pixels;
    % Calculate p_t = sum(p(i))
end

% Define the matrix A
A = [0, 0, 1, 4, 4, 5; 
     0, 1, 3, 4, 3, 4; 
     1, 3, 4, 2, 1, 3; 
     4, 4, 3, 1, 0, 0; 
     5, 4, 2, 1, 0, 0; 
     5, 5, 4, 3, 1, 0];

% Entropy-based threshold calculation
function entropy = entropy_threshold(img, threshold)
    % Calculate probabilities
    pt = calculate_pt(img);
    [~, k] = size(pt);
    eT = 0;
    % Calculate total entropy (eT)
    for i=1:k
        eT = eT + (pt(i) * log(pt(i)));
    end
    eT = -eT;
    et = 0;
    %for i=1:threshold+1
    % Calculate partial entropy up to the threshold (et)
    for i=1:threshold+1
        et = et + (pt(i) * log(pt(i)));
    end
    et = -et;
    %disp(et);
    p = 0;
    %for i=1:threshold+1
    % Calculate probability up to the threshold
    for i=1:threshold+1
        p = p + pt(i);
    end
    % Calculate entropy-based threshold
    entropy = log(p * (1-p)) + (et/p) + ((eT-et)/(1-p));
end

% Function to compute all entropy values for thresholds
function all_entropies = get_all_entropies(A)
    result = customTabulateAll(A);
    %disp(result);
    entropies = [];
    endval = max(result(:, 1));
    % Calculate entropy for each threshold
    for i=1:endval
        entropies(1,i) = entropy_threshold(A, result(i,2));
    end
    %disp(entropies);
    all_entropies = entropies;
end

% Displaying entropy-based threshold results
disp('Entropy based Thresholding');

result = customTabulateAll(A);
all_entropies = get_all_entropies(A);
disp(result(:,2)');
disp(all_entropies);
disp('----------------------------------');

% Minimum cross entropy 
% Minimum cross-entropy thresholding function
function cross_entropy = cross_entropy_threshold(A, t, pt, result)
    endval = max(result(:, 1));
    %disp(result);
    %disp(pt);
    % Initialize terms for cross-entropy calculation
    term1 = 0;
    term2 = 0;
    term3 = 0;
    
    for i=1:endval
        temp1 = 0;
        if result(i, 2) == 0
            temp1 = 0;
        else
            temp1 = (result(i, 2) * pt(i) * log10(result(i, 2)));
        end
        term1 = term1 + temp1;
    end
    % Compute term2: sum for pixels <= threshold
    for i=1:t
        total_n = 0;
        total_d = 0;
        temp2 = 0;
        for j=1:t
            total_n = total_n + (result(i,2) * result(i, 3));
            total_d = total_d + result(i, 3);
        end
        %disp(total_n);
        %disp(total_d);
        % sum(i * p(i) * log(mean(t)) , t-> threshold t = 1 
        mean = total_n/total_d;
        if mean == 0
            temp2 = 0;
        else
            temp2 = result(i, 2) * pt(i) * log10(mean);
        end
        term2 = term2 + temp2;
    end
    
    % Compute term3: sum for pixels > threshold
    total_n = 0;
    total_d = 0;
    for j=t+1:endval
            total_n = total_n + (result(j,2) * result(j, 3));
            total_d = total_d + result(j, 3);
    end
    mean2 = total_n/ total_d;
    for i=t+1:endval
        temp3 = 0;
        % sum(i * p(i) * log(mean(t)) , t-> threshold t = 1 
        if mean2 == 0
            temp3 = 0;
        else
            temp3 = result(i, 2) * pt(i) * log10(mean2);
        end
        term3 = term3 + temp3;
    end
    
    %disp(term1);
    %disp(term2);
    %disp(mean2);
    %disp(term3);
    % Final cross-entropy value
    cross_entropy = term1 - (term2 + term3);
end

% Calculating cross-entropy thresholding results
t = 1;
result = customTabulateAll(A);
%disp(result);
pt = calculate_pt(A);
cross_entropies = [];
% Computing cross-entropy for each threshold
for i=1:size(result(:,1),1)
    cross_entropy_v = cross_entropy_threshold(A, i, pt, result);
    %disp(result(i,1));
    cross_entropies = [cross_entropies, cross_entropy_v];
end

disp('Cross entropy Thresholding:');
disp(result(:,2)');
disp(cross_entropies);

% 
% Variance-based threshold calculation
function result = withinClassVariance(hist_matrix, threshold, total_values_c)
    weight_sum_bg = 0;
    mean_sum_bg = 0;
    total_values = hist_matrix(end, 1);
    %m_val_count = 0;
    % Computing background statistics
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
    
    % Compute foreground statistics
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
    % Combine weighted variances
    result = (weight_fg * variance_fg) + (weight_bg * variance_bg);
end

% Function to apply binary thresholding 
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

% function to generate gaussian filter
function gaussian_fil = gen_gaussian_filter(n, c, sigma, debug_flag)
    temp_fil = [];
    e_value = exp(1);
    k = (n-1)/2;
    for i=1:n
        for j=1:n
            temp = ((i-(k+1))^2 + (j-(k+1))^2);
            if debug_flag
                disp(i-(k+1));
                disp(j-(k+1));
                disp(temp)
            end
            temp = temp / ( 2 * sigma ^ 2);
            %temp = temp / ( 2 * sigma);
            if debug_flag
                disp(i);
                disp(j);
                disp(temp)
            end
            temp_fil(i, j) = c * e_value ^ (- temp);
        end
    end
    gaussian_fil = temp_fil;
end

% applying the gaussian filter to image
function gaussian_filtered = apply_gaussian_filter(input_img, gauss_filter, debug_flag)
    if debug_flag
        disp(input_img(1:7, 1:7));
    end
    total = size(gauss_filter, 1) * size(gauss_filter, 2);
    center = (size(gauss_filter, 1)+1)/ 2;
    [row, col] = size(gauss_filter);
    row = row - 1;
    col = col -1;
    for j=1:size(input_img, 2)-col
        for i=1:size(input_img, 1)-row
            temp = input_img(i:i+row, j:j+col);
            temp = temp.*gauss_filter;
            mean_val = round(sum(sum(temp)) / total);
            if debug_flag
                disp(temp)
                disp(mean_val)
            end
            input_img(i+center-1, j+center-1) = mean_val;
        end
    end
    gaussian_filtered = round(input_img);
end

% Function to process an image with a Gaussian filter
function gaussian_img = gen_gaussian_image(img, n, c, sigma, debug_flag)
    image_data = img; 
    image_data = double(image_data);
    if debug_flag
        display(class(image_data));
    end
    gauss_filter = gen_gaussian_filter(n, c, sigma, debug_flag);
    if debug_flag
        disp(class(gauss_filter));
        disp(gauss_filter);
    end
    gaussian_img = apply_gaussian_filter(image_data, gauss_filter, debug_flag);
end
debug_flag = 0;
% processing image

% URL of the raw image file on GitHub
url = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/sat.png';

% Specify the filename to save the downloaded file
output_filename = 'sat.png';

opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary

% Download the file with try-catch block
try
    websave(output_filename, url, opts);
    disp(['File downloaded and saved as ', output_filename]);
catch ME
    disp(['Failed to download the file. Error: ', ME.message]);
end


rgb_image = imread(output_filename);
A = rgb2gray(rgb_image);

result = customTabulateAll(A);  % Returns a matrix where each row is [value, count, percentage]
all_pixels = (result(:, 2))';
%disp(all_pixels);
variance_all = zeros(size(all_pixels));
%disp(variance_all);
total_values = size(A,1) * size(A, 2);
%disp(total_values);
%Display results
%disp(result(:, :));  % Only showing value and count columns
pixel_values_len = size(all_pixels);
for i=1:pixel_values_len(2)
    temp = withinClassVariance(result, all_pixels(i), total_values);
    variance_all(i) = temp;
end
[minValue, index] = min(variance_all);
%disp(minValue);
%disp(index);
threshold_val = all_pixels(index);
%disp(threshold_val);
img_thresholded = custThreshold(A, threshold_val);
% Displaying the outputs
subplot(2,2,1);
imshow(output_filename);
title('Original Image');
subplot(2,2,2);
imshow(img_thresholded);
title('Processed Image');

% ------------------------
%processing image 1
url1 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/surveillance.png';

% Specify the filename to save the downloaded file
output_filename1 = 'surveillance.png';

opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary
% Download the file
websave(output_filename1, url1);

disp(['File downloaded and saved as ', output_filename1]);

rgb_image1 = imread(output_filename1);
A1 = rgb2gray(rgb_image1);
%A1 = imadjust(A1);
gaussian_image1 = gen_gaussian_image(A1, 3, 1, sqrt(1000000), debug_flag);

A1 = uint8(gaussian_image1);

result1 = customTabulateAll(A1);  % Returns a matrix where each row is [value, count, percentage]
all_pixels1 = (result1(:, 2))';
%disp(all_pixels);
variance_all1 = zeros(size(all_pixels1));
%disp(variance_all);
total_values1 = size(A1,1) * size(A1, 2);
%disp(total_values);
%Display results
%disp(result(:, :));  % Only showing value and count columns
pixel_values_len1 = size(all_pixels1);
for i=1:pixel_values_len(2)
    temp = withinClassVariance(result1, all_pixels1(i), total_values1);
    variance_all1(i) = temp;
end
[minValue, index1] = min(variance_all1);
%disp(minValue);
%disp(index);
threshold_val1 = all_pixels1(index1);
%disp(threshold_val);
img_thresholded1 = custThreshold(A1, threshold_val1);
% Displaying the outputs
subplot(2,2,3);
imshow(output_filename1);
title('Original Image');
subplot(2,2,4);
imshow(img_thresholded1);
title('Processed Image');

