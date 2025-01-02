% Triangles are removed from the mesh (2D or 3D) if they have too high an aspect ratio

function filtered_mesh = Triangle_Aspect_Ratio_Filtering(triangle_mesh, point_cloud_matrix, aspect_ratio_threshold, Plot)

    %% Filtering %%
        % The vertices within the triangles        
        Triangle_Vertex_Ind_fun     = @(triangle_ind) {nchoosek(triangle_ind, 2)};
        number_triangles            = size(triangle_mesh, 1);
        triangle_vertex_ind_cell    = splitapply(Triangle_Vertex_Ind_fun, triangle_mesh', 1 : number_triangles);

        Triangle_Vertex_fun         = @(vertex_matrix) point_cloud_matrix(vertex_matrix(:, 1), :) - point_cloud_matrix(vertex_matrix(:, 2), :);
        triangle_vertex_cell        = cellfun(Triangle_Vertex_fun, triangle_vertex_ind_cell, 'UniformOutput', false);

        % The aspect ratio of each triangle
        Triangle_Vertex_Length_fun  = @(vertex_matrix) sqrt(sum(vertex_matrix.^2, 2));
        triangle_vertex_length_cell = cellfun(Triangle_Vertex_Length_fun, triangle_vertex_cell, 'UniformOutput', false);

        Triangle_Aspect_Ratio_fun   = @(triangle_vertex_lengths) prod(triangle_vertex_lengths) / (dot([1, 1, -1], triangle_vertex_lengths)*dot([-1, 1, 1], triangle_vertex_lengths)*dot([1, -1, 1], triangle_vertex_lengths));
        triangle_aspect_ratio_list  = cellfun(Triangle_Aspect_Ratio_fun, triangle_vertex_length_cell);

        % Filtering
        triangle_outliers   = triangle_aspect_ratio_list > aspect_ratio_threshold;
        filtered_mesh       = triangle_mesh(~triangle_outliers, :);

    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
            
            hold on
            grid on
            
            % The triangular meshes
            num_dim = size(point_cloud_matrix, 2);

            if num_dim == 3
                % Original
                trimesh(triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'EdgeColor', 'r', 'FaceColor', [1, 1, 1], 'FaceAlpha', 1.0, 'DisplayName', 'Original mesh');

                % Filtered
                trimesh(filtered_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'EdgeColor', 'k', 'FaceColor', [1, 1, 1], 'FaceAlpha', 1.0, 'DisplayName', 'Filtered mesh');
            
                % Axes
                xlabel('x [m]')
                ylabel('y [m]')
                zlabel('z [m]')
    
                axis equal
                view(45, 45);

            elseif num_dim == 2
                % Original
                trimesh(triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Original mesh');

                % Filtered
                trimesh(filtered_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Filtered mesh');

                % Axes
                xlabel('x [m]')
                ylabel('y [m]')
            end

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            disp('The triangular mesh has been filtered and the plot will close when a key is pressed')
            pause()

            close(1)            
        end
end
