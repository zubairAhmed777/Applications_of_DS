I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];

% I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];

S = [1, 1, 1; 1, 1, 1; 1, 1, 1];

%size(I)

for i = 1: size(I, 2) -2 % looping column
    for j = 1:size(I, 1) -2 % looping row
        %disp(I(i:i+2, j:j+2))
        max_val = max(max(I(i:i+2, j:j+2).*S)); % getting max value for dilation
        %disp(max_val)
        I(j+1, i+1) = max_val;
    end
end

%% 
disp(I);
%k = max(max(I(1:3, 1:3)));
%k
%% Dilation function
%function dilated_matrix = dilation(input_matrix, Structure_matrix)
    dilated_matrix = input_matrix;
    for i = 1: size(input_matrix, 2) -2 % looping column
        for j = 1:size(input_matrix, 1) -2 % looping row
            %disp(I(i:i+2, j:j+2))
            max_val = max(max(dilated_matrix(i:i+2, j:j+2).*Structure_matrix)); % getting max value for dilation
            %disp(max_val)
            dilated_matrix(j+1, i+1) = max_val;
        end
    end
end

I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];
S = [1, 1, 1; 1, 1, 1; 1, 1, 1];

I_dilated = dilation(I, S);
disp(I_dilated);
%% Erosion function
%function eroted_matrix = erosion(input_matrix, Structure_matrix)
    eroted_matrix = input_matrix;
    for i = 1: size(input_matrix, 2) -2 % looping column
        for j = 1:size(input_matrix, 1) -2 % looping row
            %disp(input_matrix(i:i+2, j:j+2))
            min_val = min(min(eroted_matrix(i:i+2, j:j+2).*Structure_matrix)); % getting max value for dilation
            %disp(max_val)
            eroted_matrix(i+1, j+1) = min_val;
        end
    end
end

I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];
S = [1, 1, 1; 1, 1, 1; 1, 1, 1];

I_eroted = erosion(I, S);
disp(I_eroted);

%% Opening : E --> D
%% Closing : D --> D

%% Erosion function
function eroted_matrix = erosion(input_matrix, Structure_matrix)
    eroted_matrix = input_matrix;
    for i = 1: size(input_matrix, 2) -2 % looping column
        for j = 1:size(input_matrix, 1) -2 % looping row
            %disp(input_matrix(i:i+2, j:j+2))
            min_val = min(min(eroted_matrix(i:i+2, j:j+2).*Structure_matrix)); % getting max value for dilation
            %disp(max_val)
            eroted_matrix(i+1, j+1) = min_val;
        end
    end
end

% Dilation function
function dilated_matrix = dilation(input_matrix, Structure_matrix)
    dilated_matrix = input_matrix;
    for i = 1: size(input_matrix, 2) -2 % looping column
        for j = 1:size(input_matrix, 1) -2 % looping row
            %disp(I(i:i+2, j:j+2))
            max_val = max(max(dilated_matrix(i:i+2, j:j+2).*Structure_matrix)); % getting max value for dilation
            %disp(max_val)
            dilated_matrix(j+1, i+1) = max_val;
        end
    end
end

% Opening : E --> D
function opening_matrix = opening(input_matrix, Structure_matrix)
    eroted_matrix = erosion(input_matrix, Structure_matrix);
    opening_matrix = dilation(eroted_matrix, Structure_matrix);
end

% Closing : D --> E
function closing_matrix = closing(input_matrix, Structure_matrix)
    dilated_matrix = dilation(input_matrix, Structure_matrix);
    closing_matrix = erosion(dilated_matrix, Structure_matrix);
end

%white top hat
function white_top_hat_matrix = white_top_hat(input_matrix, Structure_matrix)
    opening_mat = opening(input_matrix, Structure_matrix);
    white_top_hat_matrix = input_matrix - opening_mat;
end

%Black top hat
function black_top_hat_matrix = black_top_hat(input_matrix, Structure_matrix)
    closing_mat = closing(input_matrix, Structure_matrix);
    black_top_hat_matrix = closing_mat - input_matrix;
end

I = [10, 20, 85, 97, 55; 40,60, 70, 66, 52; 9, 70, 90, 87,12; 15, 54, 33, 60, 11; 6, 26, 73, 59, 9];
S = [1, 1, 1; 1, 1, 1; 1, 1, 1];

black_th = black_top_hat(I, S);
white_th = white_top_hat(I, S);

disp("Black Hat:")
disp(black_th);
disp("White Hat:")
disp(white_th);