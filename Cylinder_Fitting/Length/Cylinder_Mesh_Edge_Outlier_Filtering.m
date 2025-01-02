% This script removes edge points that are either on near-vertical surfaces, or not connected to others of their type    

function [point_cloud_edge_bool_filtered, axial_edge_points, neighbour_edge_outliers] = Cylinder_Mesh_Edge_Outlier_Filtering(point_cloud_edge_boolean, point_cloud_matrix, surface_triangle_mesh, cylinder_dir, alignment_threshold, edge_classes, number_edge_classes, Plot)
    
    %% Axial points %%
        % This boolean is true if a point is classified as at least one type of edge
        edge_boolean_sum    = sum(point_cloud_edge_boolean, 1);
        edge_boolean        = min(1, edge_boolean_sum);
        edge_boolean        = logical(edge_boolean);

        % Top and bottom points that are on the surface nearly parallel to the cylinder axis are removed
        % This means that the surface normal vector is approx. orthogonal
        triangle_normal_vector_matrix = Triangle_Normal_Vector(surface_triangle_mesh, point_cloud_matrix);

        % The alignment between the normal vectors and cylinder axes
        [~, ~, radial_bool, ~] = Vector_Alignment_Check(cylinder_dir, triangle_normal_vector_matrix, alignment_threshold);        
        
        % The points contained in these radial and axial triangles
        radial_points   = surface_triangle_mesh(radial_bool, :);
        radial_points   = unique(radial_points);
        
        axial_points    = surface_triangle_mesh(~radial_bool, :);
        axial_points    = unique(axial_points);

        % When the normal vectors and cylinder axis are opposed, the triangles themselves are axial
        % Only the points that are exclusively in radial triangles are desired, so exclusively axial points are considered outliers
        axial_points    = setdiff(axial_points, radial_points);
        
        % Narrowed down to edge points
        edge_points         = find(edge_boolean == true);
        axial_edge_points   = intersect(axial_points, edge_points);
        
        % Their classification is removed
        point_cloud_edge_bool_filtered                       = point_cloud_edge_boolean;
        point_cloud_edge_bool_filtered(:, axial_edge_points) = false;
        
        % This boolean is true if it is classified as at least one type of radial edge
        radial_edge_boolean_sum = sum(point_cloud_edge_bool_filtered, 1);
        radial_edge_boolean     = min(1, radial_edge_boolean_sum);
        radial_edge_boolean     = logical(radial_edge_boolean);
        radial_edge_ind         = find(radial_edge_boolean);
        num_radial_edge_points  = length(radial_edge_ind);

    %% Outliers from their neighbours %%
        neighbour_edge_outlier_bool = false(1, num_radial_edge_points);
    
        for i = 1 : num_radial_edge_points
            % The edge class(es) of this point
            ind                     = radial_edge_ind(i);
            point_edge_boolean      = point_cloud_edge_bool_filtered(:, ind);
            point_edge_classes      = find(point_edge_boolean == true);
            num_point_edge_classes  = length(point_edge_classes);
            
            % Triangles which contain this point
            [triangle_ind, ~] = find(surface_triangle_mesh == ind);
            
            % The neighbours of this point
            triangle_points_total   = surface_triangle_mesh(triangle_ind, :);
            triangle_points         = unique(triangle_points_total);
            
            neighbours = triangle_points(triangle_points ~= ind);     % The original point is not of interest
            
            % The edge classes of these points
            for c = 1 : num_point_edge_classes
                % If the edge class of the original point is not present in its neighbours, the classification is removed
                edge_class              = point_edge_classes(c);
                neighbour_edge_bool     = point_cloud_edge_bool_filtered(edge_class, neighbours);                
                frequency_edge_classes  = sum(neighbour_edge_bool, 2);
                                
                if frequency_edge_classes == 0
                    % It is marked as a neighbour outlier
                    neighbour_edge_outlier_bool(i) = true;

                    % The filtered boolean is updated
                    point_cloud_edge_bool_filtered(edge_class, ind) = false;
                end
            end
        end
        
        % The indices of the neighbour outliers
        neighbour_edge_outliers = radial_edge_ind(neighbour_edge_outlier_bool);

    %% Plot %%
        if Plot == true
            %--% Axial points and outliers %--%
            % The types of points that have been removed
            removed_points_cell     = {axial_edge_points, neighbour_edge_outliers};
            removed_points_names    = {'Axial points', 'Neighbour outliers'};
            
            % Colormap
            num_removed_types   = length(removed_points_cell);
            removed_cmap        = cbrewer('qual', 'Set1', max(num_removed_types, 3));  
            removed_cmap        = max(removed_cmap, 0);
            removed_cmap        = min(removed_cmap, 1);                                
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])
            
            for r = 1 : num_removed_types
                % The removed points of this category
                removed_points      = removed_points_cell{r};
                removed_points_name = removed_points_names{r};
                removed_colour      = removed_cmap(r, :);
                
                subplot(1, num_removed_types, r)
                hold on
                grid on

                % The triangular mesh
                trimesh(surface_triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'EdgeColor', 'k', 'FaceColor', [1, 1, 1], 'FaceAlpha', 1.0, 'DisplayName', 'Mesh');

                % The original edge points
                scatter3(point_cloud_matrix(radial_edge_boolean, 1), point_cloud_matrix(radial_edge_boolean, 2), point_cloud_matrix(radial_edge_boolean, 3), 'MarkerEdgeColor', 'b', 'DisplayName', 'Original edge points');
                
                % The removed points of this category
                scatter3(point_cloud_matrix(removed_points, 1), point_cloud_matrix(removed_points, 2), point_cloud_matrix(removed_points, 3), 'filled', 'MarkerFaceColor', removed_colour, 'DisplayName', removed_points_name);

                % Axes
                xlabel('x');
                ylabel('y');
                xlabel('z');

                axis equal

                view(45, 45)

                % Legend
                legend('show', 'location', 'westoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end
            
            %--% (Filtered) edges %--%
            % Edge colour map
            edge_cmap = cbrewer('qual', 'Set1', max(2 * number_edge_classes, 3));       % Note the factor two as each class is shown filtered and unfiltered
            edge_cmap = max(edge_cmap, 0);
            edge_cmap = min(edge_cmap, 1);
            
            figure(2)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])

            for e = 1 : number_edge_classes
                subplot(1, number_edge_classes, e);
                hold on
                grid on
            
                % The triangular mesh
                trimesh(surface_triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'EdgeColor', 'k', 'FaceColor', [1, 1, 1], 'FaceAlpha', 1.0, 'DisplayName', 'Mesh');

                % The used colours
                edge_colour             = edge_cmap(e, :);
                edge_colour_filtered    = edge_cmap(number_edge_classes + e, :);
                            
                % The identifying strings
                edge_string             = edge_classes{e};
                edge_string_filtered    = sprintf('%s %s', edge_string, '(filtered)');
                                
                % The edge points
                edge_type_boolean       = point_cloud_edge_boolean(e, :);
                scatter3(point_cloud_matrix(edge_type_boolean, 1), point_cloud_matrix(edge_type_boolean, 2), point_cloud_matrix(edge_type_boolean, 3), 'filled', 'MarkerFaceColor', edge_colour, 'DisplayName', edge_string);
            
                % The filtered edge points
                edge_boolean_filtered   = point_cloud_edge_bool_filtered(e, :);
                scatter3(point_cloud_matrix(edge_boolean_filtered, 1), point_cloud_matrix(edge_boolean_filtered, 2), point_cloud_matrix(edge_boolean_filtered, 3), 'filled', 'MarkerFaceColor', edge_colour_filtered, 'DisplayName', edge_string_filtered);
                            
                % Axes
                xlabel('x');
                ylabel('y');
                xlabel('z');

                axis equal

                view(45, 45)

                % Legend
                legend('show', 'location', 'westoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end

            % The figure is saved
            disp('The figures are saved and closed, and the script continues upon a button-press');
            pause();
            export_fig(1, 'Removed_Edge_Points.png');
            export_fig(2, 'Mesh_Edge_Points.png');
            
            close([1, 2]);            
        end
        
end