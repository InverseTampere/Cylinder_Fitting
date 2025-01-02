% This script computes the x and y coordinates of a circle, given the radius and centroid location
% Note that they start at the top of the circle, and run clockwise

function [x_circle_list, y_circle_list] = Circle_Coordinates(radius, x_centroid, y_centroid, number_coord)

    %% The coordinates %%
        % Angles at the desired discretisation level
        phi_list = linspace(0, 2*pi, number_coord);
        
        % The x and y coordinates
        x_circle_list = x_centroid + radius * sin(phi_list);
        y_circle_list = y_centroid + radius * cos(phi_list);

end