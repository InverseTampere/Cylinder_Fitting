% This script performs 2D Delaunay triangulation of a 3D cylinders

function [Delaunay_triangles_matrix, num_triangles] = Delaunay_Triangulation_Cylinder(cylinder_coords, cyl_centre, cyl_vector, cyl_radius, Triangulation_3D_Plot)

    %% Manual inputs %%
        % Filtering
        triangle_area_threshold     = 3;        % [-] When the area of a triangle is smaller or greater than the median by this factor, it is seen as an outlier
        Filter_Plot                 = false;     % [true, false] Creates a plot that shows the effect of filtering
        
        % Clustering
        Cluster_Plot                = false;    % [true, false] Shows a plot of the mesh clusters
        
        % Edge detection
        angle_margin                = 015;      % [deg] For instance, a near-vertical edge can be considered to be part of the top of the point cloud. If this is not desired,
                                                %       This margin requires a minimum tilt
        edge_margin                 = 000;      % [-]   The number of vertices of which the edge vertex may be inside of. This is recommended to be a very small fraction of
                                                %       The total number of vertices
                    
        Edge_Diagnostics            = false;    % [true, false] Creates diagnostic plots for each vertex
        Edge_Plot                   = false;    % [true, false] Creates a plot of the final result
        
    %% Point cloud unfolding %%
        % The point cloud is unfolded from 3D to 2D
        cutting_angle = 0;
        [x_unf_list, y_unf_list] = Point_Cloud_Unfolding(cylinder_coords, cyl_vector, cyl_centre, cutting_angle);
        coord_unf_matrix = [x_unf_list, y_unf_list];
        
    %% Initial Delaunay triangulation %%
        % Delaunay triangulation on unfolded point cloud
        Delaunay_triangles_matrix = delaunay(x_unf_list, y_unf_list);
    
        % Triangles that are too small or large are removed
        [Delaunay_triangles_filtered_matrix, ~, ~, ~, ~] = Mesh_Area_Filtering(Delaunay_triangles_matrix, coord_unf_matrix, triangle_area_threshold, Filter_Plot);

    %% Triangle clustering %%
        % The mesh may be disconnected after filtering, leading to multiple disconnected clusters of triangles
        [Triangle_clusters_cell, number_triangles_cluster_list, number_clusters] = Mesh_Clustering(Delaunay_triangles_filtered_matrix, coord_unf_matrix, Cluster_Plot);
    
    %% Delaunay triangulation of each cluster %%
        % The determined Delaunay triangles for each cluster with overlap
        Delaunay_triangles_cell = cell(1, number_clusters);
        number_triangles_list   = zeros(1, number_clusters);
        
        % The edges that are determined to create this overlap
        edge_point_cell         = cell(1, number_clusters);
        edge_point_top_cell     = cell(1, number_clusters);
        edge_point_right_cell   = cell(1, number_clusters);
        edge_point_bottom_cell  = cell(1, number_clusters);
        edge_point_left_cell    = cell(1, number_clusters);
        
        for c = 1 : number_clusters
            % This cluster's triangles
            Delaunay_triangles_cluster  = Triangle_clusters_cell{c};
            num_triangles               = number_triangles_cluster_list(c);            

            %--% Index translation from point cloud to cluster %--%
            cluster_points                      = unique(Delaunay_triangles_cluster);
            num_points                          = length(cluster_points);
            
            Delaunay_triangles_cluster_indices  = zeros(num_triangles, 3);
            
            for i = 1 : num_points
                % Indices that belong to this point, are given the new index
                point   = cluster_points(i);
                indices = Delaunay_triangles_cluster == point;
                
                Delaunay_triangles_cluster_indices(indices) = i;
            end
                
            % The cluster's point cloud
            cluster_coord_unf_matrix    = coord_unf_matrix(cluster_points, :);
            cluster_x_unf_list          = cluster_coord_unf_matrix(:, 1);
            cluster_y_unf_list          = cluster_coord_unf_matrix(:, 2);
            
            %--% Point cloud overlap %--%
            % As in reality a cylinder only has a top and bottom edge, the left and right edges of the unfolded plane are repeated at the opposite sides

            % The edges are detected for this cluster
            [cluster_edge_ind_list, cluster_edge_ind_top_list, cluster_edge_ind_bottom_list, cluster_edge_ind_left_list, cluster_edge_ind_right_list] = Delaunay_Triangulation_Edge_Detection(Delaunay_triangles_cluster_indices, cluster_coord_unf_matrix, angle_margin, edge_margin, Edge_Diagnostics, Edge_Plot);

            % The edges are translated according to the radius
            cluster_x_edge_list_left    = 2*cyl_radius  + cluster_x_unf_list(cluster_edge_ind_left_list);
            cluster_y_edge_list_left    = cluster_y_unf_list(cluster_edge_ind_left_list);
            cluster_x_edge_list_right   = -2*cyl_radius + cluster_x_unf_list(cluster_edge_ind_right_list);
            cluster_y_edge_list_right   = cluster_y_unf_list(cluster_edge_ind_right_list);

            % The point cloud with overlap
            cluster_x_overlap_list = [cluster_x_unf_list; cluster_x_edge_list_left; cluster_x_edge_list_right];    
            cluster_y_overlap_list = [cluster_y_unf_list; cluster_y_edge_list_left; cluster_y_edge_list_right];

            cluster_overlap_coords = [cluster_x_overlap_list, cluster_y_overlap_list];

            %--% Delaunay triangulation %--%
            % Triangles are created of this overlapping point cloud
            Delaunay_triangles_cluster_overlap_matrix = delaunay(cluster_overlap_coords);
            
            % Triangles that are too small or large are removed
            [Delaunay_triangles_cluster_overlap_matrix, ~, ~, ~, ~] = Mesh_Area_Filtering(Delaunay_triangles_cluster_overlap_matrix, cluster_overlap_coords, triangle_area_threshold, Filter_Plot);

            %--% Overlap index translation %--%
            % The indices of the overlapping segments correspond are given their cluster indices
            num_edge_points_left    = length(cluster_edge_ind_left_list);
            num_edge_points_right   = length(cluster_edge_ind_right_list);

            Delaunay_triangles_cluster_matrix = Delaunay_triangles_cluster_overlap_matrix;

            for i = 1 : num_edge_points_left
                false_index = num_points + i;
                true_index  = cluster_edge_ind_left_list(i);

                Delaunay_triangles_cluster_matrix(Delaunay_triangles_cluster_overlap_matrix == false_index) = true_index;
            end

            for j = 1 : num_edge_points_right
                false_index = num_points + num_edge_points_left + j;
                true_index  = cluster_edge_ind_right_list(j);

                Delaunay_triangles_cluster_matrix(Delaunay_triangles_cluster_overlap_matrix == false_index) = true_index;
            end

            % Duplicate triangles are removed
            Delaunay_triangles_cluster_matrix   = unique(sort(Delaunay_triangles_cluster_matrix, 2), 'rows', 'stable');
            num_triangles                       = size(Delaunay_triangles_cluster_matrix, 1);
            number_triangles_list(c)            = num_triangles;
            
            %--% Index translation from cluster to point cloud %--%
            Delaunay_triangles_cluster  = zeros(num_triangles, 3);
            
            edge_point_list             = zeros(1, length(cluster_edge_ind_list));
            edge_point_top_list         = zeros(1, length(cluster_edge_ind_top_list));
            edge_point_right_list       = zeros(1, length(cluster_edge_ind_right_list));
            edge_point_bottom_list      = zeros(1, length(cluster_edge_ind_bottom_list));
            edge_point_left_list        = zeros(1, length(cluster_edge_ind_left_list));
                        
            for i = 1 : num_points
                % Indices are given their original point back
                point = cluster_points(i);
                
                % The triangles
                triangle_indices = Delaunay_triangles_cluster_matrix == i;
                Delaunay_triangles_cluster(triangle_indices) = point;            
                
                % The edges
                edge_point_list(cluster_edge_ind_list == i)                 = point;
                edge_point_top_list(cluster_edge_ind_top_list == i)         = point;
                edge_point_right_list(cluster_edge_ind_right_list == i)     = point;
                edge_point_bottom_list(cluster_edge_ind_bottom_list == i)   = point;
                edge_point_left_list(cluster_edge_ind_left_list == i)       = point;
            end
            
            % The triangle matrix is appended
            Delaunay_triangles_cell{c} = Delaunay_triangles_cluster;
            
            % The edges are appended
            edge_point_cell{c}          = edge_point_list;
            edge_point_top_cell{c}      = edge_point_top_list;
            edge_point_right_cell{c}    = edge_point_right_list;
            edge_point_bottom_cell{c}   = edge_point_bottom_list;
            edge_point_left_cell{c}     = edge_point_left_list;
        end
        
    %% Combined results %%
        % The Delaunay triangles of each cluster
        Delaunay_triangles_matrix   = vertcat(Delaunay_triangles_cell{:});
        num_triangles               = size(Delaunay_triangles_matrix, 1);
        
        % The edges used to creat the overlap
        edge_point_ind          = horzcat(edge_point_cell{:});
        edge_point_top_ind      = horzcat(edge_point_top_cell{:});
        edge_point_right_ind    = horzcat(edge_point_right_cell{:});
        edge_point_bottom_ind   = horzcat(edge_point_bottom_cell{:});
        edge_point_left_ind     = horzcat(edge_point_left_cell{:});
        
    %% Plot %%
    if Triangulation_3D_Plot == true
        % Boolean to indicate edges or not
        num_points          = length(x_unf_list);
        edge_point_boolean  = false(1, num_points);
        edge_point_boolean(edge_point_ind) = true;
        
        % Edge colour map
        edge_cmap = cbrewer('qual', 'Set1', 4);
        edge_cmap = max(edge_cmap, 0);
        edge_cmap = min(edge_cmap, 1);
        
        figure(1)
        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
        set(gcf, 'color', [1, 1, 1])

        hold on
        grid on

        % The mesh
        trimesh(Delaunay_triangles_matrix, cylinder_coords(:, 1), cylinder_coords(:, 2), cylinder_coords(:, 3), 'EdgeColor', 'b', 'DisplayName', 'Delaunay mesh');
        
        % The point cloud
        scatter3(cylinder_coords(~edge_point_boolean, 1), cylinder_coords(~edge_point_boolean, 2), cylinder_coords(~edge_point_boolean, 3), 'filled', 'm', 'DisplayName', 'Point cloud');

        % The edge points
        scatter3(cylinder_coords(edge_point_left_ind, 1), cylinder_coords(edge_point_left_ind, 2), cylinder_coords(edge_point_left_ind, 3), 'filled', 'MarkerFaceColor', edge_cmap(4, :), 'DisplayName', 'Left edge');
        scatter3(cylinder_coords(edge_point_right_ind, 1), cylinder_coords(edge_point_right_ind, 2), cylinder_coords(edge_point_right_ind, 3), 'filled', 'MarkerFaceColor', edge_cmap(2, :), 'DisplayName', 'Right edge');
        scatter3(cylinder_coords(edge_point_top_ind, 1), cylinder_coords(edge_point_top_ind, 2), cylinder_coords(edge_point_top_ind, 3), 'filled', 'MarkerFaceColor', edge_cmap(1, :), 'DisplayName', 'Top edge');
        scatter3(cylinder_coords(edge_point_bottom_ind, 1), cylinder_coords(edge_point_bottom_ind, 2), cylinder_coords(edge_point_bottom_ind, 3), 'filled', 'MarkerFaceColor', edge_cmap(3, :), 'DisplayName', 'Bottom edge');

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
        
        % Saving the figure
        disp('The figure is saved and closed when a button is pressed, and the script continues.')
        pause()
        
        export_fig(1, 'Delaunay_Triangulation_3D.png');
        close(1);
    end    
end