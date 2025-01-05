% Defining a 5x5 matrix I (input image)
I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];

% I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];
 
% Define a 3x3 averaging filter S
S = [1, 1, 1;
    1, 1, 1;
    1, 1, 1];

% Mean Filter
% This function applies a mean filter to the input image
function mean_fil = mean_filter_(input_img, input_filter)
    for i=1:size(input_img, 1)-2
        for j=1:size(input_img, 2) -2
            temp = input_img(i:i+2, j:j+2).*input_filter;
            %disp(temp)
            s = sum(sum(temp));
            %disp(s);
            input_img(i+1, j+1) = round(s/sum(input_filter(:)));
        end
    end
    mean_fil = input_img;
end
%disp(I);
% Display the output of the Mean Filter
disp('Mean Filter Output:')
disp(mean_filter_(I, S));

% Median Filter
% This function applies a median filter to the input image
function median_fil = median_filter(input_img, input_filter)
    for i=1:size(input_img, 1) -2
        for j=1:size(input_img, 2) -2
            %sorted = sort(input_img(i:i+2, j:j+2));
            %disp(sorted);
            %disp(input_img(i:i+2, j:j+2));
            conv_mat = input_img(i:i+2, j:j+2) .* input_filter;
            vector=[];
            for z=1:size(conv_mat,1)
                vector = [vector conv_mat(z,:)];
            end
            vector = sort(vector);
            median_ele = vector(5);
            input_img(i+1, j+1) = median_ele;
        end
    end
    median_fil = input_img;
end
% Displaying the output of the Median Filter
disp('Median Filter Output:')
disp(median_filter(I, S));
% Gaussian filter

% sigma sq = 3
% n = 3, c = 1
% g[i, j] = c * e ^ -(i^2 + j^2)/ 2 * sigma ^ 2
% function to generate gaussian filter
function gaussian_fil = gen_gaussian_filter(n, c, sigma, debug_flag)
    temp_fil = [];
    e_value = exp(1); % Euler's number
    k = (n-1)/2; % Calculate center offset
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
            temp_fil(i, j) = c * e_value ^ (- temp); % Calculate Gaussian value
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
    image_data = imread(img); % Read the input image
    image_data = double(image_data);
    if debug_flag
        display(class(image_data));
    end
    gauss_filter = gen_gaussian_filter(n, c, sigma, debug_flag); % Generate Gaussian filter
    if debug_flag
        disp(class(gauss_filter));
        disp(gauss_filter);
    end
    gaussian_img = apply_gaussian_filter(image_data, gauss_filter, debug_flag); % Apply the Gaussian filter
end
debug_flag = 0;

disp('*********** Start **********');
disp(' ');

url = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/noise_1.jpg';

% Specify the filename to save the downloaded file
image_filename1 = 'noise_1.jpg';

% Download and process the first image with try-catch block
try
    websave(image_filename1, url);
    disp(['File downloaded and saved as ', image_filename1]);
catch ME
    disp(['Failed to download the file. Error: ', ME.message]);
end

gaussian_image1 = gen_gaussian_image(image_filename1, 3, 1, sqrt(1000000), debug_flag);

gaussian_image1 = uint8(gaussian_image1);

subplot(2, 2, 1);
imshow(image_filename1);
title(['Original image ', image_filename1]);
subplot(2, 2, 2);
imshow(gaussian_image1);
title(['Gaussian filter to ', image_filename1]);

% URL of the raw image file on GitHub
url2 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/noise_2.jpg';

% Specify the filename to save the downloaded file
image_filename2 = 'noise_2.jpg';

% Download and process the second image with try-catch block
try
    websave(image_filename2, url2);
    disp(['File downloaded and saved as ', image_filename2]);
catch ME
    disp(['Failed to download the file. Error: ', ME.message]);
end


gaussian_image2 = gen_gaussian_image(image_filename2, 3, 1, sqrt(1000000), debug_flag);

gaussian_image2 = uint8(gaussian_image2);
subplot(2, 2, 3);
imshow(image_filename2);
title(['Original image ', image_filename2]);

subplot(2, 2, 4);
imshow(gaussian_image2);
title(['Gaussian filter to ', image_filename2]);
%imshow(image_1);
disp('*********** Done **********');

