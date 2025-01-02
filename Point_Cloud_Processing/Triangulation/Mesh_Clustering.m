% This script clusters a triangular mesh, and returns the clusters in a cell array
% Note that coordinates (in 2D) need only be given if plots are desired

function [Triangle_clusters_cell, number_triangles_list, number_clusters] = Mesh_Clustering(Triangle_matrix, coordinate_matrix, Plot)
    
    %% Loop to create clusters %%
        % This loop iteratively clusters the triangles until no triangles remain
        Triangle_matrix_original = Triangle_matrix;
        
        Triangle_matrix = sort(Triangle_matrix, 2);
        num_triangles   = size(Triangle_matrix, 1);
        
        Triangle_clusters_cell  = cell(1, num_triangles);
        number_triangles_list   = zeros(1, num_triangles);
        
        number_clusters = 0;
        while ~isempty(Triangle_matrix)
            number_clusters = number_clusters + 1;
            
            % A triangle is selected at random from the remaining ones
            num_triangles   = size(Triangle_matrix, 1);
            random_triangle = randi([1, num_triangles], 1);
            centre_triangle = Triangle_matrix(random_triangle, :);
            
            % It is used to initiate a new cluster
            % Triangles that are connected to this cluster are found iteratively
            Triangle_cluster            = centre_triangle;
            Triangle_cluster_prior      = Triangle_cluster;
            
            Cluster_Complete = false;
            
            while Cluster_Complete == false
                % Neighbouring triangles are found for each triangle                
                num_cluster_triangles_prior = size(Triangle_cluster_prior, 1);
                Triangle_cluster_post       = cell(1, num_cluster_triangles_prior);
                
                for t = 1 : num_cluster_triangles_prior
                    % This triangle
                    triangle = Triangle_cluster_prior(t, :);
                    
                    % Neighbours are defined as triangles which have one or more points in common
                    neighbour_cell = cell(1, 3);

                    for i = 1 : 3
                        % Triangles which have this point in common
                        point_index = triangle(i);

                        [triangle_ind, ~]   = find(Triangle_matrix == point_index);
                        neighbour_triangles = Triangle_matrix(triangle_ind, :);
                        neighbour_cell{i}   = neighbour_triangles;
                    end
                    
                    % All neighbours are assigned to the cluster
                    neighbours = vertcat(neighbour_cell{:});
                    
                    Triangle_cluster_post{t} = neighbours;
                end

                % Triangles that are new
                Triangle_cluster_post   = vertcat(Triangle_cluster_post{:});
                Triangle_cluster_post   = sort(Triangle_cluster_post, 2);
                
                Triangle_cluster_new    = setdiff(Triangle_cluster_post, Triangle_cluster, 'rows', 'stable');
                Triangle_cluster_new    = sort(Triangle_cluster_new, 2);
                
                % The cluster thus far
                Triangle_cluster        = union(Triangle_cluster_new, Triangle_cluster, 'rows', 'stable');
                Triangle_cluster        = sort(Triangle_cluster, 2);
                
                % If the cluster has not increased in size, it is complete
                if isempty(Triangle_cluster_new)
                    % Plot showing this cluster
                    if Plot == true
                        figure(number_clusters)
                        % Set the size and white background color
                        set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                        set(gcf, 'color', [1, 1, 1])

                        hold on
                        grid on
                        
                        % The triangle meshes
                        triplot(Triangle_matrix_original, coordinate_matrix(:, 1), coordinate_matrix(:, 2), 'color', 'b', 'DisplayName', 'Full mesh');
                        triplot(Triangle_cluster, coordinate_matrix(:, 1), coordinate_matrix(:, 2), 'color', 'g', 'DisplayName', 'Cluster mesh');

                        % The point cloud
                        scatter(coordinate_matrix(:, 1), coordinate_matrix(:, 2), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Point cloud');
                        
                        % Axes
                        xlabel('x');
                        ylabel('y');

                        axis equal

                        % Legend
                        legend('show', 'location', 'eastoutside');

                        set(gca, 'FontSize', 10);
                        set(gca, 'LineWidth', 2);

                        hold off
                    end
                    
                    Cluster_Complete                            = true;
                    Triangle_clusters_cell{number_clusters}     = Triangle_cluster;
                    
                    num_triangles_cluster                       = size(Triangle_cluster, 1);
                    number_triangles_list(number_clusters)      = num_triangles_cluster;
                    
                    % The triangles appended to this cluster are removed
                    [~, ind_cluster_triangles, ~]               = intersect(Triangle_matrix, Triangle_cluster, 'rows', 'stable');
                    Triangle_matrix(ind_cluster_triangles, :)   = [];
                    
                % Otherwise, the cluster is updated and another search will take place only over the newly found triangles
                else
                    Triangle_cluster_prior = Triangle_cluster_new;                    
                end
            end
        end
        
        % Empty clusters are removed
        empty_ind                           = cellfun(@isempty, Triangle_clusters_cell);
        
        Triangle_clusters_cell(empty_ind)   = [];
        number_triangles_list(empty_ind)    = [];
        
        % The figures are closed
        if Plot == true
            disp('The figures are closed and the script continues upon a button-press');
            pause();

            close(1 : number_clusters);
        end
end