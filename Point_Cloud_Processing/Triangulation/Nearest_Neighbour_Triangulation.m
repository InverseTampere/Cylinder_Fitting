% This script performs k-nearest neighbour 2D triangulation of the given 3D point cloud
% The number of neighbours is based off of the total number of points and the given neighbour percentage

function Triangle_matrix = Nearest_Neighbour_Triangulation(coordinate_matrix, neighbour_percentage, Triangulation_Diagnostics, Triangulation_Plot)

    %% Triangulation %%
        % The number of neighbours 
        num_points      = size(coordinate_matrix, 1);
        num_neighbours  = round(num_points * neighbour_percentage / 100);
        
        num_neighbours  = min(num_points - 1, num_neighbours);
        
        if Triangulation_Diagnostics == true
            fprintf('For kNN triangulation, %g neighbours are used \n', num_neighbours);
        end

        % The triangles with k nearest neighbours are found separately for each point
        Triangle_cell = cell(1, num_points);

        for i = 1 : num_points
             % This point's coordinates
            point = coordinate_matrix(i, :);

            % Its nearest neighbours
            offset_matrix = coordinate_matrix - point;
            distance_list = sqrt(sum(offset_matrix.^2, 2));

            [~, neighbour_order] = sort(distance_list, 'ascend');

            neighbours = neighbour_order(2 : num_neighbours + 1);       % Note that the first neighbour is the point itself, and is thus ignored

            % The triangles formed by the combinations of neighbours
            neighbour_combinations_total    = nchoosek(neighbours, 2);
            neighbour_combinations          = unique(sort(neighbour_combinations_total, 2), 'rows', 'stable');

            number_combinations = size(neighbour_combinations, 1);
            index_list          = repmat(i, [number_combinations, 1]);
            triangles           = [index_list, neighbour_combinations];     % The triangles include the original point index
            Triangle_cell{i}    = triangles;       
            
            % Plot showing the neighbours and the generated triangles
            if Triangulation_Diagnostics == true
                figure(1)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])

                hold on
                grid on            

                % The generated mesh for this point
                trimesh(triangles, coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'EdgeColor', 'b', 'DisplayName', 'kNN mesh');

                % The point cloud
                scatter3(coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'filled', 'm', 'DisplayName', 'Point cloud');

                % The central point
                scatter3(point(1), point(2), point(3), 'filled', 'c', 'DisplayName', 'Central point');
                
                % The neighbours
                scatter3(coordinate_matrix(neighbours, 1), coordinate_matrix(neighbours, 2), coordinate_matrix(neighbours, 3), 'filled', 'r', 'DisplayName', 'Neighbours');
                                
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
                disp('The figure is closed when a button is pressed, and the script continues.');
                pause();
                close(1);
            end
        end
                
        % The total set of triangles
        Triangle_matrix = vertcat(Triangle_cell{:});
        Triangle_matrix = unique(sort(Triangle_matrix, 2), 'rows', 'stable');
        
    %% Plot %%
        if Triangulation_Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])

            hold on
            grid on            
            
            % The generated mesh
            trimesh(Triangle_matrix, coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'EdgeColor', 'b', 'DisplayName', 'kNN mesh');
            
            % The point cloud
            scatter3(coordinate_matrix(:, 1), coordinate_matrix(:, 2), coordinate_matrix(:, 3), 'filled', 'm', 'DisplayName', 'Point cloud');
            
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

            export_fig(1, 'kNN_Triangulation.png');
            close(1);
        end
end