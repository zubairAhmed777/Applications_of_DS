% Function to pad an image with a 1-pixel border
function result = custPad(img)
    [row, col] = size(img);
    result = ones(row+2, col+2);
    [r_row, r_col] = size(result);
    for i=2:r_row-1
        for j=2:r_col-1
            result(i,j) = img(i-1, j-1);
        end
    end
end

% function to do sum
function total = customSum(array)
    % Initialize total to 0
    total = 0;
    % Loop through each element of the array
    for i = 1:length(array)
        total = total + array(i);
    end
end

% Function to perform horizontal convolution
function result = hxConv(temp, fil)
    % Extract middle row and apply filter
    filteredRow = temp(2,:) .* fil;
    % Use custom sum function
    result = customSum(filteredRow);
end

% Function to perform vertical convolution
function result = hyConv(temp, fil)
    % Extract middle column and apply filter
    filteredCol = temp(:,2) .* fil;
    % Use custom sum function
    result = customSum(filteredCol);
end

% Function to compute gradients (magnitude and angle)
function [mag, ang] = cust_hog(A, pA, fx, fy)
    mag = zeros(size(A)); % Initialize magnitude array
    ang = zeros(size(A)); % Initialize angle array
    [row, col] = size(pA); % Size of padded image
    for i=1:row-2
        for j=1:col-2
            temp = pA(i:i+2, j:j+2); % Extract 3x3 region
            hx = hxConv(temp, fx); % Horizontal gradient
            hy = hyConv(temp, fy); % Vertical gradient
            gi = sqrt((hx * hx) + (hy * hy)); % Gradient magnitude
            y_rad = atan2(hy, hx); % Gradient angle in radians
            ai = rad2deg(y_rad); % Convert angle to degrees
            if ai < 0
                ai = 180 + ai; % Normalize angle to [0, 180]
            end
            mag(i,j) = gi; % Store magnitude
            ang(i,j) = floor(ai); % Store angle
            %{
            if isnan(ai)
                disp(hx);
                disp(hy);
            end
            %}
        end
    end
end

% Function to generate histogram of oriented gradients (HoG)
function [index1, weight1, index2, weight2] = genHOG(angle, magnitude)
    numBins = 9;
    % Generate bin edges between 0 and 180 degrees
    binEdges = linspace(0, 180, numBins + 1);

    % Determine which bin the value belongs to
    binIndex = discretize(angle, binEdges);

    % Default lower and upper bounds
    lowerBound = -1;
    upperBound = -1;

    if ~isnan(binIndex)
        % Get the lower and upper bounds of the bin
        lowerBound = binEdges(binIndex);
        upperBound = binEdges(binIndex + 1);
    end

    % Calculate weights for the two bins
    weight1 = round((abs(angle - upperBound) / (180 / numBins)) * magnitude);
    weight2 = round((abs(angle - lowerBound) / (180 / numBins)) * magnitude);

    % Assign indices
    index1 = binIndex;

    % Handle wrap-around for the next bin index
    if binIndex + 1 > numBins
        index2 = 1; % Wrap around to the first bin
    else
        index2 = binIndex + 1;
    end
end

% Function to generate HoG features for an image block
function hog_temp = genHOG_features(magnitudes, angles)
    % Initialize HoG temp array
    hog_temp = zeros(1, 9);
    
    [rows, columns] = size(angles);
    for i=1:rows
        for j=1:columns
            % Input values
            angle = angles(i,j);
            magnitude = magnitudes(i,j);
            
            % Call genHOG function
            [index1, weight1, index2, weight2] = genHOG(angle, magnitude);
            
            
            % Update the HoG feature vector
            hog_temp(1, index1) = hog_temp(1, index1) + weight1;
            hog_temp(1, index2) = hog_temp(1, index2) + weight2;
        end
    end
end

A = [0,0,0,1,4,4,5,5;
    0,0,0,1,4,4,5,5; 
    0,0,1,3,4,3,4,4; 
    1,1,3,4,2,1,3,3; 
    4,4,4,3,1,0,0,0; 
    5,5,4,2,1,0,0,0; 
    5,5,5,4,3,1,0,0;
    5,5,5,4,3,1,0,0];

fx = [-1,0,1];
fy = fx';
pA = custPad(A);
%disp(pA);
% Compute magnitudes and angles
[magnitudes, angles] = cust_hog(A, pA, fx, fy);
%{
disp('Magnitude')
disp(magnitudes);
disp('Angles');
disp(angles);
%}
% Generate HoG features for the sample matrix
hog_temp = genHOG_features(magnitudes, angles);
disp('Image A:')
disp(A);
disp('Hog Features for matrix A :');
disp(hog_temp);

% Process template image for HoG features
url1 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/template.png';
% Specify the filename to save the downloaded file
output_filename1 = 'template.png';
opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjusting timeout

% Download the file with try-catch block
try
    websave(output_filename1, url1, opts);
    disp(['Template File downloaded and saved as --> ', output_filename1]);
catch ME
    disp(['Failed to download the template file. Error: ', ME.message]);
end

% Read and preprocess template image
% Read the RGB image
rgbImage = imread(output_filename1); 

% Convert the RGB image to grayscale
grayImage = rgb2gray(rgbImage);
%resizing the image to 128 * 64
resizedImage = imresize(grayImage, [128, 64]);
% Get the size of the image
[rows, cols] = size(resizedImage);

fx = [-1,0,1];
fy = fx';

all_hog_features = [];

% Dividing the image into 16x16 blocks and compute HoG features
for i = 1:8:(rows-8)
    temp_hog_features = [];
    for j = 1:8:(cols-8)
        % Calculate the block indices
        temp = resizedImage(i:i+15, j:j+15);
        %disp([i, i+15, " ---", j, j+15] );
        % splitting the block into 4 cells of 8 * 8 size
        %gen hog for 1st cell
        temp1 = temp(1:8,1:8);
        ptemp1 = custPad(temp1);
        [magnitudes1, angles1] = cust_hog(temp1, ptemp1, fx, fy);
        hog_temp1 = genHOG_features(magnitudes1, angles1);
        % gen hog for 2nd cell
        temp2 = temp(1:8,9:16);
        ptemp2 = custPad(temp2);
        [magnitudes2, angles2] = cust_hog(temp2, ptemp2, fx, fy);
        hog_temp2 = genHOG_features(magnitudes2, angles2);
        %gen hog for 3rd cell
        temp3 = temp(9:16,1:8);
        ptemp3 = custPad(temp3);
        [magnitudes3, angles3] = cust_hog(temp3, ptemp3, fx, fy);
        hog_temp3 = genHOG_features(magnitudes3, angles3);
        %gen hog for 4th cell
        temp4 = temp(9:16,9:16);
        ptemp4 = custPad(temp4);
        [magnitudes4, angles4] = cust_hog(temp4, ptemp4, fx, fy);
        hog_temp4 = genHOG_features(magnitudes4, angles4);
        temp_hog_features = [temp_hog_features, hog_temp1, hog_temp2, hog_temp3, hog_temp4];
    end
    %disp('----------------');
    all_hog_features = [all_hog_features; temp_hog_features];
end

%disp(all_hog_features);
disp('HoG Features are generated for image & saved in "all_hog_features" variable');




