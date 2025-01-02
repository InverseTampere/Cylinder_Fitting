% This script determines the points belonging to the top and bottom edges of a cylinder's point cloud using a triangular mesh
% Whether or not a point belongs to the side or the top/bottom is determined via its triangles' alignments with the cylinder direction
% A low threshold is strict. It is expected to be given in degrees

function [edge_classes, number_edge_classes, point_cloud_edge_bool_f, point_cloud_edge_cell_f, number_edge_points, surface_triangle_mesh] = Cylinder_Mesh_Edge_Points(alignment_threshold, cylinder_dir, cylinder_centre, Point_Cloud_Distributions, Diagnostics, Plot)

    %% Inputs %%
        % Distributions
        num_points              = Point_Cloud_Distributions.number_distributions;
        distribution_mu_cell    = Point_Cloud_Distributions.distribution_mu_cell;

    %% Mesh creation %%
        % Delaunay triangulation to create a mesh of 3D triangles
        point_cloud_matrix                  = vertcat(distribution_mu_cell{:});
        [surface_triangle_mesh, num_facets] = Surface_Mesh_Creator_DT(point_cloud_matrix);

    %% Intersections of the mesh parallel to the cylinder axis %%
        % The cylinder axis vector placed at each point
        cylinder_vector_matrix  = repmat(cylinder_dir, [num_points, 1]);
        
        % The intersections with the edge facets
        delta_matrix = NaN(num_points, num_facets);
        
        for f = 1 : num_facets
            % The coordinates of this triangle
            triangle_ind            = surface_triangle_mesh(f, :);
            triangle_coord_matrix   = point_cloud_matrix(triangle_ind, :);
            
            % The intersects of points that are part of this triangle are not of interest
            point_cloud_matrix_f                    = point_cloud_matrix;
            point_cloud_matrix_f(triangle_ind, :)   = NaN;
            
            % The distances to the intersections between the axial vectors and this triangle
            [~, delta_list, ~, ~]   = Triangle_Vector_Intersection(triangle_coord_matrix, point_cloud_matrix_f, cylinder_vector_matrix);
            delta_matrix(:, f)      = delta_list;
        end
        
    %% Edge classification %%
        % Edges are detected and classified in the order Top - Bottom
        edge_classes        = {'Top', 'Bottom'};
        number_edge_classes = length(edge_classes);
    
        point_cloud_edge_boolean = false(number_edge_classes, num_points);
        
        % Edges have intersections in either zero or one direction
        delta_sign_matrix       = sign(delta_matrix);
        num_post_intersects     = sum(delta_sign_matrix > 0, 2);
        num_prior_intersects    = sum(delta_sign_matrix < 0, 2);
        
        % Additionally top and bottom points must lie above and below the cylinder centre respectively
        [~, height_list, ~] = Point_to_Vector_Projection(point_cloud_matrix, cylinder_dir, cylinder_centre);
    
        % Top points
        top_points = num_post_intersects == 0 & height_list > 0;
        point_cloud_edge_boolean(1, top_points) = true;
        
        % Bottom points
        bottom_points = num_prior_intersects == 0 & height_list < 0;
        point_cloud_edge_boolean(2, bottom_points)  = true;
        
        % Plot for diagnostic purposes
        if Diagnostics == true
            % Edge colour map
            edge_cmap = cbrewer('qual', 'Set1', max(number_edge_classes, 3));
            edge_cmap = max(edge_cmap, 0);
            edge_cmap = min(edge_cmap, 1);
            
            for i = 1 : num_points
                % Only points with intersects are shown
                delta_list  = delta_matrix(i, :);
                delta_list  = delta_list(~isnan(delta_list));
    
                if isempty(delta_list)
                    continue
                end
    
                % Message to show the type of edge
                disp('-------------');
                edge_boolean = point_cloud_edge_boolean(:, i);
    
                if ~edge_boolean
                    disp('This point is not considered an edge');
                else
                    edge_string = edge_classes(edge_boolean);
                    edge_string = strjoin(edge_string, '/');
    
                    fprintf('This point is considered a %s edge \n', edge_string);
                end
    
                % The intersects
                point               = point_cloud_matrix(i, :);
                intersections       = point + delta_list' * cylinder_dir;
                intersections_post  = intersections(delta_list > 0, :);
                intersections_prior = intersections(delta_list <= 0, :);
    
                % The point's colour
                if ~isempty(find(edge_boolean == true, 1))
                    point_colour = edge_cmap(edge_boolean, :);
                    point_colour = mean(point_colour, 1);
                else
                    point_colour = 'b';
                end
    
                figure(i)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])                
    
                hold on
                grid on
    
                % The surface facets
                trisurf(surface_triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'LineWidth', 2, 'FaceColor', [1, 1, 1], 'FaceAlpha', 0.5, 'DisplayName', 'Bounding facets');
    
                % The point
                scatter3(point(1), point(2), point(3), 'filled', 'MarkerFaceColor', point_colour, 'DisplayName', 'Point');
    
                % The axis direction
                scaled_cyl_dir = cyl_length * cyl_dir / norm(cyl_dir);
                plot3(point(1) + scaled_cyl_dir(1) * [-1, 1], point(2) + scaled_cyl_dir(2) * [-1, 1], point(3) + scaled_cyl_dir(3) * [-1, 1], 'LineWidth', 2, 'color', point_colour, 'DisplayName', 'Axis vector');
                
                % The intersections
                scatter3(intersections_post(:, 1), intersections_post(:, 2), intersections_post(:, 3), 'filled', 'MarkerFaceColor', 'c', 'DisplayName', 'Post intersects');
                scatter3(intersections_prior(:, 1), intersections_prior(:, 2), intersections_prior(:, 3), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'Prior intersects');
    
                % Axes
                xlabel('x');
                ylabel('y');
                zlabel('z');
    
                axis equal
    
                % Viewing angle
                view(45, 45);
    
                % Legend
                legend('show', 'location', 'eastoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off    
    
                % Pause
                disp('The script will continue and the plot will close when a button is pressed.');
                pause();
                close(i);
            end
        end

    %% Outlier filtering %%
        % Outlier edge points are removed
        [point_cloud_edge_bool_f, ~, ~] = Cylinder_Mesh_Edge_Outlier_Filtering(point_cloud_edge_boolean, point_cloud_matrix, surface_triangle_mesh, cylinder_dir, alignment_threshold, edge_classes, number_edge_classes, Plot);
        
        point_cloud_edge_cell_f = cell(1, number_edge_classes);

        for e = 1 : number_edge_classes
            edge_bool                   = point_cloud_edge_bool_f(:, e);
            point_cloud_edge_matrix     = point_cloud_matrix(edge_bool, :);
            point_cloud_edge_cell_f{e}  = point_cloud_edge_matrix;
        end

        number_edge_points = sum(point_cloud_edge_bool_f, 2);

end