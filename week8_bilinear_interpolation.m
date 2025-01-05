% Original Image
originalImage = [10, 20; 30, 40];

% Original grid points
originalX = [1, 2];
originalY = [1, 2];

% Defining finer grid points (4x4 output)
fineX = linspace(1, 2, 4);
fineY = linspace(1, 2, 4);

% Initializing the output matrix
interpolatedOutput = zeros(length(fineY), length(fineX));

% Performing bilinear interpolation manually
for rowIndex = 1:length(fineY)
    for colIndex = 1:length(fineX)
        % Defining the coordinates for interpolation
        interpX = fineX(colIndex);
        interpY = fineY(rowIndex);

        % Identifying the surrounding points
        xLower = 1; xUpper = 2;
        yLower = 1; yUpper = 2;

        % Values at the four surrounding points
        valueQ11 = originalImage(1, 1);
        valueQ12 = originalImage(2, 1);
        valueQ21 = originalImage(1, 2);
        valueQ22 = originalImage(2, 2);

        % Bilinear interpolation formula
        tempValue = valueQ11 * ((xUpper - interpX) * (yUpper - interpY)) + ...
                    valueQ21 * ((interpX - xLower) * (yUpper - interpY)) + ...
                    valueQ12 * ((xUpper - interpX) * (interpY - yLower)) + ...
                    valueQ22 * ((interpX - xLower) * (interpY - yLower));
        interpolatedOutput(rowIndex, colIndex) = floor(tempValue);
    end
end

% Displaying results
disp('Original Matrix:');
disp(originalImage);
disp(['Size of Original Matrix: ', num2str(size(originalImage))]);

disp('Resulting Interpolated Matrix (4x4):');
disp(interpolatedOutput);
disp(['Size of Interpolated Matrix: ', num2str(size(interpolatedOutput))]);

%
% Bilinear Interpolation Function
function resizedImage = bilinear_interpolation(inputImage, scaleFactor)
    % Getting the original dimensions
    [numRows, numCols] = size(inputImage);

    % Calculating the new dimensions
    newNumRows = floor(numRows * scaleFactor);
    newNumCols = floor(numCols * scaleFactor);

    % Initialize the resized image
    resizedImage = zeros(newNumRows, newNumCols);

    % Scaling factors for coordinates
    rowScale = numRows / newNumRows;
    colScale = numCols / newNumCols;

    % Iterate over each pixel in the new image
    for newRowIndex = 1:newNumRows
        for newColIndex = 1:newNumCols
            % Finding the corresponding position in the original image
            originalX = newColIndex * colScale;
            originalY = newRowIndex * rowScale;

            % Identifying the surrounding points
            xLower = floor(originalX); xUpper = ceil(originalX);
            yLower = floor(originalY); yUpper = ceil(originalY);

            % Handling edge cases
            xLower = max(1, xLower);
            xUpper = min(numCols, xUpper);
            yLower = max(1, yLower);
            yUpper = min(numRows, yUpper);

            % Pixel values at the surrounding points
            valueQ11 = double(inputImage(yLower, xLower));
            valueQ12 = double(inputImage(yUpper, xLower));
            valueQ21 = double(inputImage(yLower, xUpper));
            valueQ22 = double(inputImage(yUpper, xUpper));

            % Interpolation weights
            xWeight = originalX - xLower;
            yWeight = originalY - yLower;

            % Bilinear interpolation
            interpolatedValue = ...
                valueQ11 * (1 - xWeight) * (1 - yWeight) + ...
                valueQ21 * xWeight * (1 - yWeight) + ...
                valueQ12 * (1 - xWeight) * yWeight + ...
                valueQ22 * xWeight * yWeight;

            % Assigning the interpolated value to the new image
            resizedImage(newRowIndex, newColIndex) = interpolatedValue;
        end
    end

    % Converting to uint8 for display
    resizedImage = uint8(resizedImage);
end

% Grayscale Conversion Function
function grayscaleImage = convert_to_grayscale(inputImage)
    % Check if the image is RGB
    if ndims(inputImage) == 3 && size(inputImage, 3) == 3
        [numRows, numCols, ~] = size(inputImage);
        grayscaleImage = zeros(numRows, numCols);

        % Convert each pixel to grayscale
        for rowIndex = 1:numRows
            for colIndex = 1:numCols
                redChannel = inputImage(rowIndex, colIndex, 1);
                greenChannel = inputImage(rowIndex, colIndex, 2);
                blueChannel = inputImage(rowIndex, colIndex, 3);

                % Calculate grayscale intensity
                grayscaleImage(rowIndex, colIndex) = 0.2989 * redChannel + ...
                                                     0.5870 * greenChannel + ...
                                                     0.1140 * blueChannel;
            end
        end

        grayscaleImage = uint8(grayscaleImage);
        disp('Image was RGB and has been converted to Grayscale.');
    else
        grayscaleImage = inputImage;
        disp('Image is already Grayscale or single-channel. No conversion needed.');
    end
end

% Example Usage - Image 1
imageURL = 'https://github.com/zubairAhmed777/Applications_of_DS/raw/main/colour.jpg';
try
    websave('colour.jpg', imageURL);
    disp('Image successfully downloaded.');
catch
    error('Failed to fetch the image. Check the URL or internet connection.');
end

originalColourImage = imread('colour.jpg');
figure;
imshow(originalColourImage);
title('Original Image');

grayscaleImage1 = convert_to_grayscale(originalColourImage);
figure;
imshow(grayscaleImage1);
title('Grayscale Converted Image');

disp(['Size of Original Image: ', num2str(size(grayscaleImage1))]);

scaleFactor = 2;
resizedImage1 = bilinear_interpolation(grayscaleImage1, scaleFactor);
disp(['Size of Resized Image: ', num2str(size(resizedImage1))]);

figure;
imshow(resizedImage1);
title('Resized Image (Bilinear Interpolation)');

% Example Usage - Image 2
% Define the URL of the image
imageURL = 'https://github.com/zubairAhmed777/Applications_of_DS/raw/main/gray.jpg';
try
    websave('gray.jpg', imageURL);
    disp('Image successfully downloaded.');
catch
    error('Failed to fetch the image. Check the URL or internet connection.');
end

originalGrayImage = imread('gray.jpg');
figure;
imshow(originalGrayImage);
title('Original Image');

grayscaleImage2 = convert_to_grayscale(originalGrayImage);
figure;
imshow(grayscaleImage2);
title('Grayscale Converted Image');

disp(['Size of Original Image: ', num2str(size(grayscaleImage2))]);

resizedImage2 = bilinear_interpolation(grayscaleImage2, scaleFactor);
disp(['Size of Resized Image: ', num2str(size(resizedImage2))]);

figure;
imshow(resizedImage2);
title('Resized Image (Bilinear Interpolation)');