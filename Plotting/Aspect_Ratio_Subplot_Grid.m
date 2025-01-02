% The optimal number of rows and columns of subplots is determined here for the given number of subplots and aspect ratio (w/h)
% If the aspect ratio is not given, it is assumed to be 16/9

function [number_rows, number_columns] = Aspect_Ratio_Subplot_Grid(number_subplots, aspect_ratio)

    % If the aspect ratio is not given, it is defined as 16/9
    if nargin == 1
        aspect_ratio = 16/9;
    end

    % The number of columns and rows of subplots
    number_columns  = ceil(sqrt(aspect_ratio * number_subplots));       
    number_rows     = floor(1/aspect_ratio * number_columns);       
end