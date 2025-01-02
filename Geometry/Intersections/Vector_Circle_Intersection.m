% This script computes the intersection between the given vectors and the circle described by its centroid and radius
% The vectors start at [x_start_list, y_start_list] and the direction is given by [vector_x_list, vector_y_list]

% The closest intersections are given in a separate matrix, with a boolean denoting the vectors that intersected
% All intersects, whether complex or real, are given in the full intersects matrix with their corresponding distances in the distance matrix

function [closest_intersects_matrix, full_intersects_matrix, distance_matrix, intersection_bool] = Vector_Circle_Intersection(points_matrix, vector_matrix, circle_centre, circle_radius, Diagnostics)

    %% The intersection points %%
        % The vectors are normalised
        vector_matrix = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));
        
        % The distances to the intersection points
        number_vectors          = size(points_matrix, 1);
        circle_centre_matrix    = repmat(circle_centre, [number_vectors, 1]);       % The centres are repeated for the dot product
        
        B_list = 2*dot(points_matrix - circle_centre_matrix, vector_matrix, 2);
        C_list = sum(points_matrix.^2, 2) + sum(circle_centre_matrix.^2, 2) - circle_radius^2 - 2*dot(points_matrix, circle_centre_matrix, 2);
        
        delta_matrix = 1/2 * (-B_list + [-1, 1] .* sqrt(B_list.^2 - 4*C_list));
        
        % Resulting intersection points
        full_intersects_matrix  = points_matrix + cat(3, delta_matrix(:, 1) .* vector_matrix, delta_matrix(:, 2) .* vector_matrix);
        
        % The nearest deltas
        distance_matrix = abs(delta_matrix);
        closest_ind     = (distance_matrix == min(distance_matrix, [], 2))';        % The transpose is used to simplify the next set of code
        delta_matrix    = delta_matrix';
        
        % If they are equidistant, the first one is selected
        closest_ind_sum     = sum(closest_ind, 1);
        equidistants        = closest_ind_sum == 2;
        num_equidistants    = sum(equidistants);
        
        closest_ind(:, equidistants) = repmat([1; 0], 1, num_equidistants);
        
        % Checking they are real
        intersection_bool = imag(delta_matrix(closest_ind)) == 0;        

        % The resulting nearest, real intersections
        closest_intersects_matrix   = points_matrix + delta_matrix(closest_ind) .* vector_matrix;
        closest_intersects_matrix   = closest_intersects_matrix(intersection_bool, :);
                
    %% Diagnostics plot %%
    if Diagnostics == true
        % The number of coordinates used for the circle
        number_coord = 1e3;

        figure(1)
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
        set(gcf, 'color', [1, 1, 1])   
        
        hold on
        grid on

        % The circle
        [x_circle_list, y_circle_list] = Circle_Coordinates(circle_radius, circle_centre(1), circle_centre(2), number_coord);
        plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Circle');

        % The points
        scatter(points_matrix(:, 1), points_matrix(:, 2), 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Points');

        % The vectors
        scale = max([real(distance_matrix(:)); 2*circle_radius], [], 'all');

        for v = 1 : number_vectors
            point           = points_matrix(v, :);
            vector          = vector_matrix(v, :);
            scaled_vector   = scale * vector;

            pl_vec = plot(point(1) + [0, scaled_vector(1)], point(2) + [0, scaled_vector(2)], 'LineWidth', 1, 'color', 'r', 'DisplayName', 'Vectors');

            if v > 1
                pl_vec.HandleVisibility = 'Off';
            end
        end

        % All intersects
        for i = 1 :2
            intersects_matrix = squeeze(full_intersects_matrix(intersection_bool, :, i));
            scatter(intersects_matrix(:, 1), intersects_matrix(:, 2), 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Intersects');
        end

        % The closest intersects
        scatter(closest_intersects_matrix(:, 1), closest_intersects_matrix(:, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'DisplayName', 'Closest intersects');

        % Axes
        axis equal

        xlabel('x [m]');
        ylabel('y [m]');

        % Legend
        legend('show', 'location', 'northoutside');
        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 2);

        % Pause message
        disp('The script will continue and the figure will close when a button is pressed');
        pause();
        close(1);
    end    
end