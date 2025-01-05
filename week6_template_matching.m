% Define an example input matrix and template
I = [0 0 1 4 4 5;
0 1 3 4 3 4;
1 3 4 2 1 3;
4 4 3 1 0 0;
5 4 2 1 0 0;
5 5 4 3 1 0];

T = [1 0 0;
1 0 0;
3 1 0];

% Download template image
url1 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/template.png';
% Specifying the filename to save the downloaded file
output_filename1 = 'template.png';
opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary
% Downloading the file with try-catch block
try
    websave(output_filename1, url1, opts);
    disp(['Template File downloaded and saved as ', output_filename1]);
catch ME
    disp(['Failed to download the template file. Error: ', ME.message]);
end

% Download target image
url2 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/snap_4.png';
% Specifying the filename to save the downloaded file
output_filename2 = 'snap_4.png';
opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary
% Downloading the file with try-catch block
try
    websave(output_filename2, url2, opts);
    disp(['Image File downloaded and saved as ', output_filename2]);
catch ME
    disp(['Failed to download the image file. Error: ', ME.message]);
end

% Download video file
url3 = 'https://raw.githubusercontent.com/zubairAhmed777/Applications_of_DS/main/video.mp4';
% Specifying the filename to save the downloaded file
output_filename3 = 'video.mp4';
opts = weboptions('Timeout', 15, 'CertificateFilename', ''); % Adjust timeout if necessary
% Downloading the file with try-catch block
try
    websave(output_filename3, url3, opts);
    disp(['Video File downloaded and saved as ', output_filename3]);
catch ME
    disp(['Failed to download the video file. Error: ', ME.message]);
end

% Read the images
image = imread(output_filename2);
template = imread(output_filename1);

% Convert template to grayscale if it is RGB
if size(template, 3) == 3
        template = rgb2gray(template);
end

[tmplRows, tmplCols] = size(template);
% Template matching using SAD
disp('Processing for SAD');
[bestMatchRow, bestMatchCol, sadMap] = templateMatchingSAD(image, template);
%fprintf('Best Match Location: Row = %d, Col = %d\n', bestMatchRow, bestMatchCol);

%Displaying the image and bounding boxes with SAD 
figure;
imshow(uint8(image)), title('Template Matching Result Using SAD:');
hold on;
rectangle('Position', [bestMatchCol, bestMatchRow, tmplCols, tmplRows], ...
              'EdgeColor', 'b', 'LineWidth', 2);
hold off;
% Template matching using SSD
disp('Processing for SSD');
[bestMatchRow2, bestMatchCol2, ssdMap] = templateMatchingSSD(image, template);
figure;
imshow(uint8(image)), title('Template Matching Result Using SSD:');
hold on;
rectangle('Position', [bestMatchCol2, bestMatchRow2, tmplCols, tmplRows], ...
              'EdgeColor', 'y', 'LineWidth', 2);
hold off;
% Template matching using NCC
disp('Processing for NCC');
[bestMatchRow3, bestMatchCol3, nccMap] = templateMatchingNCC(image, template);
figure;
imshow(uint8(image)), title('Template Matching Result Using NCC');
hold on;
rectangle('Position', [bestMatchCol3, bestMatchRow3, tmplCols, tmplRows], ...
              'EdgeColor', 'r', 'LineWidth', 2);
hold off;

% Process video file with template matching and save the result
%templateMatchingVideo('video.mp4', 'template.png', 'video_ncc.mp4');

% Template Matching using SAD
function [bestMatchRow, bestMatchCol, sadMap] = templateMatchingSAD(image, template)
    % Convert input to grayscale if they are RGB
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    
    if size(template, 3) == 3
        template = rgb2gray(template);
    end
    
    % Convert to double for computation
    image = double(image);
    template = double(template);
    
    % Get dimensions of the image and the template
    [imgRows, imgCols] = size(image);
    [tmplRows, tmplCols] = size(template);
    
    % Initialize SAD map to store SAD values
    sadMap = inf(imgRows - tmplRows + 1, imgCols - tmplCols + 1);
    
    % Slide the template over the image
    for k = 1:(imgRows - tmplRows + 1)
        for l = 1:(imgCols - tmplCols + 1)
            % Compute SAD using explicit formula
            sadValue = 0;
            for i = 1:tmplRows
                for j = 1:tmplCols
                    % Compute absolute difference
                    diff = abs(template(i, j) - image(i + k - 1, j + l - 1));
                    sadValue = sadValue + diff;
                end
            end
            
            % Store the SAD value in the map
            sadMap(k, l) = sadValue;
        end
    end
    
    % Find the position of the minimum SAD value
    [~, minIdx] = min(sadMap(:));
    [bestMatchRow, bestMatchCol] = ind2sub(size(sadMap), minIdx);
    
    % Display results
    %fprintf('Best Match Location: Row = %d, Col = %d\n min val = %d\n', bestMatchRow, bestMatchCol, minIdx);
    
    % Visualize the match on the original image
    %{
    figure;
    imshow(uint8(image)), title('Template Matching Result');
    hold on;
    rectangle('Position', [bestMatchCol, bestMatchRow, tmplCols, tmplRows], ...
              'EdgeColor', 'y', 'LineWidth', 2);
    hold off;
    %}
end

% Template Matching using SSD
function [bestMatchRow, bestMatchCol, ssdMap] = templateMatchingSSD(image, template)
    % Convert input to grayscale if they are RGB
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    
    if size(template, 3) == 3
        template = rgb2gray(template);
    end
    
    % Convert to double for computation
    image = double(image);
    template = double(template);
    
    % Get dimensions of the image and the template
    [imgRows, imgCols] = size(image);
    [tmplRows, tmplCols] = size(template);
    
    % Initialize SSD map to store SSD values
    ssdMap = inf(imgRows - tmplRows + 1, imgCols - tmplCols + 1);
    
    % Slide the template over the image
    for k = 1:(imgRows - tmplRows + 1)
        for l = 1:(imgCols - tmplCols + 1)
            % Compute SSD using explicit formula
            ssdValue = 0;
            for i = 1:tmplRows
                for j = 1:tmplCols
                    % Compute squared difference
                    diff = (template(i, j) - image(i + k - 1, j + l - 1))^2;
                    ssdValue = ssdValue + diff;
                end
            end
            
            % Store the SSD value in the map
            ssdMap(k, l) = ssdValue;
        end
    end
    
    % Find the position of the minimum SSD value
    [~, minIdx] = min(ssdMap(:));
    [bestMatchRow, bestMatchCol] = ind2sub(size(ssdMap), minIdx);
    %{
    % Display results
    fprintf('Best Match Location: Row = %d, Col = %d\n', bestMatchRow, bestMatchCol);
    
    % Visualize the match on the original image
    figure;
    imshow(uint8(image)), title('Template Matching Result');
    hold on;
    rectangle('Position', [bestMatchCol, bestMatchRow, tmplCols, tmplRows], ...
              'EdgeColor', 'r', 'LineWidth', 2);
    hold off;
    %}
end

% Template Matching using NCC
function [bestMatchRow, bestMatchCol, nccMap] = templateMatchingNCC(image, template)
    % Convert input to grayscale if they are RGB
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    
    if size(template, 3) == 3
        template = rgb2gray(template);
    end
    
    % Convert to double for computation
    image = double(image);
    template = double(template);
    
    % Get dimensions of the image and the template
    [imgRows, imgCols] = size(image);
    [tmplRows, tmplCols] = size(template);
    
    % Precompute the squared template sum for normalization
    templateSum = sum(template(:));
    templateSumSq = sum(template(:).^2);
    
    % Initialize NCC map to store NCC values
    nccMap = -inf(imgRows - tmplRows + 1, imgCols - tmplCols + 1);
    
    % Slide the template over the image
    for k = 1:(imgRows - tmplRows + 1)
        for l = 1:(imgCols - tmplCols + 1)
            % Extract region of interest (ROI) from the image
            roi = image(k:(k + tmplRows - 1), l:(l + tmplCols - 1));
            
            % Compute the numerator (cross-correlation)
            numerator = sum(sum(template .* roi));
            
            % Compute the denominator (normalization factor)
            roiSumSq = sum(roi(:).^2);
            denominator = sqrt(templateSumSq * roiSumSq);
            
            % Avoid division by zero
            if denominator ~= 0
                nccMap(k, l) = numerator / denominator;
            else
                nccMap(k, l) = 0; % Assign a very low value
            end
        end
    end
    
    % Find the position of the maximum NCC value
    [~, maxIdx] = max(nccMap(:));
    [bestMatchRow, bestMatchCol] = ind2sub(size(nccMap), maxIdx);
    
    % Display results
    %fprintf('Best Match Location: Row = %d, Col = %d\n', bestMatchRow, bestMatchCol);
    
    % Visualize the match on the original image
    %{
    figure;
    imshow(uint8(image)), title('Template Matching Result');
    hold on;
    rectangle('Position', [bestMatchCol, bestMatchRow, tmplCols, tmplRows], ...
              'EdgeColor', 'r', 'LineWidth', 2);
    hold off;
    %}
end

% Template Matching in Video
function templateMatchingVideo(videoFile, templateFile, outputVideoFile)
    % Read the input video file
    disp(['Video file is : ', videoFile]);
    disp(['Template file is : ', templateFile]);
    disp(['Output Video file is : ', outputVideoFile]);
    
    video = VideoReader(videoFile);
    
    % Read the template image
    template = imread(templateFile);
    if size(template, 3) == 3
        template = rgb2gray(template);
    end
    template = double(template);
    
    % Preparing output video writer
    outputVideo = VideoWriter(outputVideoFile, 'MPEG-4');
    open(outputVideo);
    
    % Processing each frame
    while hasFrame(video)
        % Read a frame
        frame = readFrame(video);
        frameGray = rgb2gray(frame);
        frameGray = double(frameGray);
        
        % Apply NCC for template matching
        [bestRow, bestCol, nccMap] = templateMatchingNCC(frameGray, template);
        
        % Annotate the matched region on the frame
        annotatedFrame = insertShape(uint8(frame), 'Rectangle', ...
            [bestCol, bestRow, size(template, 2), size(template, 1)], ...
            'Color', 'green', 'LineWidth', 3);
        
        % Write annotated frame to the output video
        writeVideo(outputVideo, annotatedFrame);
        
        % Display the processed frame
        imshow(annotatedFrame);
        title('Template Matching Result Video');
        drawnow;
    end
    
    % Closing the video writer
    close(outputVideo);
    disp(['Processed video saved to: ', outputVideoFile]);    
end