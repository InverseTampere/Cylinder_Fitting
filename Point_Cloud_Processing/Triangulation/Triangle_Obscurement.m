% This script returns an unobscured triangle mesh w.r.t. the scanner of interest
% Obscurement is tested through intersecting vectors from the scanner

function [surface_triangle_matrix, num_surface_triangles, obscured_triangles_boolean] = Triangle_Obscurement(bounding_triangles, coordinate_matrix, scanner_loc, Plot)

    %% Triangle information %%
        % Triangle centroids
        num_facets = size(bounding_triangles, 1);
        
        triangle_centroid_func  = @(triangle_ind) {mean(coordinate_matrix(triangle_ind, :), 1)};
        centroid_cell           = splitapply(triangle_centroid_func, bounding_triangles', 1 : num_facets);    
        centroid_matrix         = vertcat(centroid_cell{:});

        % Vectors from the scanner to these centroids
        scanner_loc_matrix          = repmat(scanner_loc, [num_facets, 1]);
        scanner_triangle_vectors    = centroid_matrix - scanner_loc;
        scanner_triangle_vectors    = scanner_triangle_vectors ./ sqrt(sum(scanner_triangle_vectors.^2, 2));
    
    %% Intersections %%
        % Distances for all vectors to intersect each bounding triangle
        delta_matrix = NaN(num_facets);

        for f = 1 : num_facets
            % This triangle's coordinates
            triangle_ind            = bounding_triangles(f, :);
            triangle_coord_matrix   = coordinate_matrix(triangle_ind, :);

            % The distances to the intersections
            [~, delta_list, ~, ~]   = Triangle_Vector_Intersection(triangle_coord_matrix, scanner_loc_matrix, scanner_triangle_vectors);
            delta_matrix(f, :)      = delta_list;            
        end

        % The minimum distance for each vector to the triangle intersections should be for their own triangle
        min_delta_list  = min(delta_matrix, [], 1);
        diag_delta_list = diag(delta_matrix)';

        obscured_triangles_boolean = min_delta_list ~= diag_delta_list;

        % Obscured triangles are removed
        surface_triangle_matrix = bounding_triangles(~obscured_triangles_boolean, :);
        num_surface_triangles   = size(surface_triangle_matrix, 1);

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])
        
            hold on
            grid on

            % The bounding facets
            trimesh(bounding_triangles, coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'EdgeColor', 'k', 'FaceAlpha', 0.5, 'DisplayName', 'Bounding triangles');

            % The surface triangles
            trimesh(surface_triangle_matrix, coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Surface triangles');
            
            % Point cloud
            scatter3(coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'DisplayName', 'Point cloud');

            % Axes
            xlabel('x');
            ylabel('y');
            zlabel('z');

            axis equal

            % Viewing angle
            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            % Pausing message
            disp('The triangle obscurement script will end and the plot will close upon a button-press');
            pause();
            close(1);
        end
end