%% 
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
%% 
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
    
    if isnan(variance_bg)
        variance_bg = 0;
    end

    if isnan(variance_fg)
        variance_fg = 0;
    end


    variance_fg = variance_fg/weight_sum_fg;
    disp('weight bg');
    disp(weight_bg);
    disp('var bg: ')
    disp(variance_bg);
    disp('weight fg');
    disp(weight_fg)
    disp('var fg');
    disp(variance_fg);
    result = (weight_fg * variance_fg) + (weight_bg * variance_bg);
end

%% 
A = [0,0,1,4,4,5; 0,1,3,4,3,4; 1,3,4,2,1,3; 4,4,3,1,0,0; 5,4,2,1,0,0; 5,5,4,3,1,0];
result = customTabulateAll(A);  % Returns a matrix where each row is [value, count, percentage]

total_values = size(A,1) * size(A, 2);
% Display results
%disp(result(:, :));  % Only showing value and count columns
disp(withinClassVariance(result, 0, total_values));


