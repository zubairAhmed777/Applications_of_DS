I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];

% I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];

S = [1, 1, 1;
    1, 1, 1;
    1, 1, 1];

function mean_fil = mean_filter(input_img, input_filter)
    for i=1:size(input_img, 1)-2
        for j=1:size(input_img, 2) -2
            temp = input_img(i:i+2, j:j+2).*input_filter;
            %disp(temp)
            s = sum(sum(temp));
            %disp(s);
            input_img(i+1, j+1) = round(s/9);
        end
    end
    mean_fil = input_img;
end
%disp(I);
disp(mean_filter(I, S));
%% 

function median_fil = median_filter(input_img, input_filter)
    for i=1:size(input_img, 1) -2
        for j=1:size(input_img, 2) -2
            %sorted = sort(input_img(i:i+2, j:j+2));
            %disp(sorted);
            %disp(input_img(i:i+2, j:j+2));
            conv_mat = input_img(i:i+2, j:j+2) .* input_filter;
            vector=[];
            for z=1:size(conv_mat,1)
                vector=[vector conv_mat(z,:)];
            end
            median_ele = vector(5);
            %disp(median_ele);
            input_img(i+1, j+1) = median_ele;
        end
    end
    median_fil = input_img;
end

disp(median_filter(I, S));

%% Gaussian filter

% sigma sq = 3
% n = 3, c = 1
% g[i, j] = c * e ^ -(i^2 + j^2)/ 2 * sigma ^ 2

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

function gaussian_img = gen_gaussian_image(img, n, c, sigma, debug_flag)
    image_data = imread(img); 
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

disp('*********** Start **********');
disp(' ');


gaussian_image1 = gen_gaussian_image("noise_1.jpg", 3, 1, sqrt(1000000), debug_flag);

gaussian_image1 = uint8(gaussian_image1);
imshow(gaussian_image1);

if debug_flag
    image_1 = imread("noise_1.jpg");
    disp(gaussian_image1(955:960, 1275:1280));
    disp(image_1(955:960, 1275:1280));
    disp(size(gaussian_image1));
    disp(size(image_1))
end
%imshow(image_1);
disp('*********** Done **********');

%% 
image_1 = imread("noise_1.jpg");
disp(size(image_1))

disp(image_1(1:5, 1:5))


