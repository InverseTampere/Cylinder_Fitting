% This script computes the intersection point between a triangle and the given vectors
% For vectors parallel to the triangle, or intersections outside the triangle, the intersection and delta is defined as NaN

function [intersection_matrix, delta_list, parallel_vectors, outside_triangle_intersections] = Triangle_Vector_Intersection(triangle_coord_matrix, vector_start_matrix, vector_matrix)

    %% Manual inputs %%
        Plot    = false;     % [true, false] Shows the triangle, vectors and their intersections

    %% Plane on which the triangle lies %%
        % The three corner points
        point_a = triangle_coord_matrix(1, :);
        point_b = triangle_coord_matrix(2, :);
        point_c = triangle_coord_matrix(3, :);
        
        % The vector normal to the triangle
        triangle_normal_vector = cross(point_b - point_a, point_c - point_a);
        
    %% Intersections of the vectors with this plane %%
        [intersection_matrix, delta_list, parallel_vectors] = Plane_Vector_Intersection(point_a, triangle_normal_vector, vector_start_matrix, vector_matrix);
    
    %% Check if intersections lie within the triangles %%
        % They lie within the triangle if the sum of the three triangles the intersection points form with the original triangle points is equal to the area of the original triangle
        
        % Vectors from the intersections to the triangle points
        vector_u_matrix = point_a - intersection_matrix;
        vector_v_matrix = point_b - intersection_matrix;
        vector_w_matrix = point_c - intersection_matrix;
        
        % Areas of the intersections with the triangle points
        matrix_norm             = @(x) sqrt(sum(x.^2, 2));
        A_intersections_list    = 1/2 * (matrix_norm(cross(vector_u_matrix, vector_v_matrix, 2)) + matrix_norm(cross(vector_u_matrix, vector_w_matrix, 2)) + matrix_norm(cross(vector_v_matrix, vector_w_matrix, 2)));
        
        % Area of the triangle
        A_triangle      = 1/2 * norm(triangle_normal_vector);        

        % Check if they are equal
        rounding_margin = 1e-6 * A_triangle;    % As they will never be exactly equal
        A_diff_list     = abs(A_intersections_list - A_triangle);
        
        outside_triangle_intersections = A_diff_list > rounding_margin;
        
        % Intersections outside triangles are marked as NaN
        intersection_matrix(outside_triangle_intersections, :)  = NaN;
        delta_list(outside_triangle_intersections)              = NaN;
        
    %% Plot %%
        if Plot == true
            % Triangle loop
            triangle_loop = [triangle_coord_matrix; triangle_coord_matrix(1, :)];
            
            % Triangle coordinate amplitude
            UB_list = max(triangle_coord_matrix, [], 1, 'omitnan');
            LB_list = min(triangle_coord_matrix, [], 1, 'omitnan');
            
            lambda  = sqrt(sum((UB_list - LB_list).^2, 2)) / 1e1;

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
            
            hold on
            grid on
            
            % The triangle
            plot3(triangle_loop(:, 1), triangle_loop(:, 2), triangle_loop(:, 3), 'b', 'DisplayName', 'Triangle');
            
            % The normal vector to the triangle
            triangle_centroid = mean(triangle_coord_matrix, 1);
            plot3(triangle_centroid(1) + [0, triangle_normal_vector(1)], triangle_centroid(2) + [0, triangle_normal_vector(2)], triangle_centroid(3) + [0, triangle_normal_vector(3)], 'color', 'b', 'LineWidth', 2, 'DisplayName', 'Triangle normal');
    
            % The intersection points
            scatter3(intersection_matrix(:, 1), intersection_matrix(:, 2), intersection_matrix(:, 3), 'filled', 'MarkerFaceColor' ,'r', 'DisplayName', 'Intersections');
            
            % The vectors
            num_vectors = size(vector_matrix, 1);
            
            for v = 1 : num_vectors
                % The data of this vector
                vector          = vector_matrix(v, :); 
                vector_start    = vector_start_matrix(v, :);
                delta           = delta_list(v);
                
                if isnan(delta)
                    delta = lambda;
                end
                    
                pl = plot3(vector_start(1) + delta*vector(1)*[0, 1], vector_start(2) + delta*vector(2)*[0, 1], vector_start(3) + delta*vector(3)*[0, 1], 'color', 'r', 'LineWidth', 1, 'DisplayName', 'Vectors');
                
                if v > 1
                    pl.HandleVisibility = 'Off';
                end
            end
            
            % Axes
            xlim([LB_list(1), UB_list(1)] + lambda*[-1, 1]);
            ylim([LB_list(2), UB_list(2)] + lambda*[-1, 1]);
            zlim([LB_list(3), UB_list(3)] + lambda*[-1, 1]);

            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            view(45, 45);
            
            % Legend
            legend('show', 'location', 'eastoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Time to look at the figure
            disp('The figure is closed when a button is pressed, and the script continues.');
            pause();
            close(1);
        end
end
    