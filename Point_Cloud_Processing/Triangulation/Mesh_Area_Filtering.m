% This script filters out triangles that are too small or large
% Note that this means that the resulting mesh may not be fully connected, in which case the separate clusters are separate entries in the output cell array

function [Triangles_filtered_matrix, num_triangles, too_small_triangles, too_large_triangles, triangle_outliers] = Mesh_Area_Filtering(Triangles_matrix, coordinate_matrix, triangle_area_threshold, Plot)

    %% The area of each triangle %%
        num_triangles = size(Triangles_matrix, 1);
        
        triangle_area_list = zeros(1, num_triangles);
        
        for t = 1 : num_triangles
            % The coordinates of this triangle
            triangle_ind    = Triangles_matrix(t, :);
            triangle_coord  = coordinate_matrix(triangle_ind, :);
            
            point_i = triangle_coord(1, :);
            point_j = triangle_coord(2, :);
            point_k = triangle_coord(3, :);
            
            % The vertex between points i and j
            vector_t = point_i - point_j;
            
            % The projection of point k onto this vector
            proj_k = point_j + dot(point_k - point_j, vector_t) / norm(vector_t)^2 * vector_t;
            
            % The vector that passes from the projected point to k, and is normal to vector_t
            vector_n = point_k - proj_k;
            
            % The area from the projected point to points i and j
            area_proj_i = 1/2 * norm(vector_n) * norm(point_i - proj_k);
            area_proj_j = 1/2 * norm(vector_n) * norm(proj_k - point_j);
            
            % The resulting area depends on whether or not the projected point lies outside the triangle
            if      norm(proj_k - point_i) > norm(point_j - point_i) && norm(proj_k - point_i) > norm(proj_k - point_j)         % The projected point lies beyond point j
                triangle_area = area_proj_i - area_proj_j;
            elseif  norm(proj_k - point_j) > norm(point_j - point_i) && norm(proj_k - point_j) > norm(proj_k - point_i)         % The projected point lies before point i
                triangle_area = area_proj_j - area_proj_i;
            else
                triangle_area = area_proj_i + area_proj_j;
            end            
            
            triangle_area_list(t) = triangle_area;
        end
        
   %% Outlier identification and removal %%
        % Too small or large triangles are identified as outliers
        triangle_area_list      = triangle_area_list / mean(triangle_area_list);        % Normalised
        
        median_triangle_area    = median(triangle_area_list);
        triangle_area_ratios    = triangle_area_list / median_triangle_area;
        
        too_large_triangles     = triangle_area_ratios > triangle_area_threshold;
        too_small_triangles     = triangle_area_ratios < 1/triangle_area_threshold;
        
        triangle_outliers       = logical(too_large_triangles + too_small_triangles);

        % The outlier triangles are removed
        Triangles_filtered_matrix  = Triangles_matrix(~triangle_outliers, :);
        num_triangles                       = size(Triangles_filtered_matrix, 1);
        
    %% Plot %%
        if Plot == true            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])           

            hold on
            grid on
            
            % The final triangles
            triplot(Triangles_filtered_matrix, coordinate_matrix(:, 1), coordinate_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Resulting mesh');
            
            % The filtered triangles
            trisurf(Triangles_matrix(triangle_outliers, :), coordinate_matrix(:, 1), coordinate_matrix(:, 2), zeros(size(coordinate_matrix, 1), 1), 'FaceAlpha', 0.5, 'FaceColor', 'b', 'EdgeColor', 'b', 'LineWidth', 2, 'DisplayName', 'Filtered triangles');

            % The original point cloud
            scatter(coordinate_matrix(:, 1), coordinate_matrix(:, 2), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'Point cloud');

            % Axes
            xlabel('x');
            ylabel('y');
            
            axis equal
            
            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
    
            % Saving the figure
            disp('The figure is saved and closed when a button is pressed, and the script continues.')
            pause()

            export_fig(1, 'Area_Based_Triangulation_Filtering.png');
            close(1);
        end
        
end