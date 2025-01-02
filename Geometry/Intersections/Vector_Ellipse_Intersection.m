% This script computes the intersection point between a vector and an ellipse

% The closest intersections are given in a separate matrix, with a boolean denoting the vectors that intersected
% All intersects, whether complex or real, are given in the full intersects matrix with their corresponding distances in the distance matrix

function [closest_intersects_matrix, full_intersects_matrix, distance_matrix, intersection_bool] = Vector_Ellipse_Intersection(points_matrix, vector_matrix, ellipse_centre, ellipse_radii, ellipse_axes, Diagnostics)
    
    %% Transformation to a circle %%
        % The vectors are normalised
        vector_matrix   = vector_matrix ./ sqrt(sum(vector_matrix.^2, 2));
    
        % The ellipse centre is used to centre the coordinate frame
        point_matrix_c  = points_matrix - ellipse_centre;
    
        % The frame is rotated s.t. the ellipse axes align with the coordinate frame, and then the ellipse is squeezed to a unit circle
        num_dim                 = length(ellipse_centre);
        squeeze_matrix          = eye(num_dim) ./ ellipse_radii;
        
        point_matrix_r  = (squeeze_matrix * ellipse_axes * point_matrix_c')';
        vector_matrix_r = (squeeze_matrix * ellipse_axes * vector_matrix')';
        
    %% Intersections with circle %%
        unit_circle_radius  = 1;
        circle_centre       = zeros(1, num_dim);
        [closest_circle_intersects_matrix, full_circle_intersects_matrix, circle_distance_matrix, intersection_bool] = Vector_Circle_Intersection(point_matrix_r, vector_matrix_r, circle_centre, unit_circle_radius, Diagnostics);
        
    %% Transformation to ellipse intersections %%
        % Transformation back to the original coordinate frame
        expand_matrix               = ellipse_radii .* eye(num_dim);
        closest_intersects_matrix_c = (ellipse_axes' * expand_matrix * closest_circle_intersects_matrix')';
        full_intersects_matrix_c    = cat(3, (ellipse_axes' * expand_matrix * full_circle_intersects_matrix(:, :, 1)')', (ellipse_axes' * expand_matrix * full_circle_intersects_matrix(:, :, 2)')');
        
        % Centering at the origin
        closest_intersects_matrix   = closest_intersects_matrix_c + ellipse_centre;
        full_intersects_matrix      = full_intersects_matrix_c + ellipse_centre;
        
        % The distances are multiplied by the squeezing factor
        squeezing_factors   = sqrt(sum(vector_matrix_r.^2, 2));
        distance_matrix     = circle_distance_matrix ./ squeezing_factors;
        
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

        % The ellipse
        ellipse_coord_matrix = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, number_coord);
        plot(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Ellipse');

        % The vectors
        num_vectors = size(points_matrix, 1);
        scale       = max([distance_matrix(intersection_bool, :); ellipse_radii], [], 'all');

        for v = 1 : num_vectors
            point           = points_matrix(v, :);
            vector          = vector_matrix(v, :);
            scaled_vector   = scale * vector;

            pl_vec = plot(point(1) + [0, scaled_vector(1)], point(2) + [0, scaled_vector(2)], 'LineWidth', 1, 'color', 'r', 'DisplayName', 'Vectors');

            if v > 1
                pl_vec.HandleVisibility = 'Off';
            end
        end

        % All intersects
        for i = 1 : 2
            intersects_matrix = squeeze(full_intersects_matrix(intersection_bool, :, i));
            sc_int = scatter(intersects_matrix(:, 1), intersects_matrix(:, 2), 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Intersects');

            if i == 2
                sc_int.HandleVisibility = 'Off';
            end
        end

        % The closest intersects
        scatter(closest_intersects_matrix(intersection_bool, 1), closest_intersects_matrix(intersection_bool, 2), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', 'DisplayName', 'Closest intersects');

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