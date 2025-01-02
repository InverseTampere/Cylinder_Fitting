% The given matrix is projected onto the ellipse using the shortest paths
% Note that this script is faster than the analytical version

function [projected_points_matrix, distance_list] = Point_to_Ellipse_Projection_Numerical(ellipse_centre, ellipse_radii, ellipse_axes, points_matrix, number_samples, Plot)

    %% Transformation to ellipse coordinate frame %%
        % The points are placed in the ellipse-centered and oriented frame
        points_matrix_t     = points_matrix - ellipse_centre;
        points_matrix_r     = (ellipse_axes * points_matrix_t')';

        % Additionally, they are placed in the first quarter to reduce the number of samples needed
        sign_matrix         = sign(points_matrix_r);
        points_matrix_q     = sign_matrix .* points_matrix_r;
    
    %% Projection %%
        % Samples in the first quarter of the ellipse
        theta_list          = linspace(0, pi/2, number_samples)';
        ellipse_samples     = ellipse_radii .* [cos(theta_list), sin(theta_list)];
        
        % Distance from each point to the samples
        Distance_fun        = @(x, y) sqrt(sum((ellipse_samples - [x, y]).^2, 2));
        distance_cell       = arrayfun(Distance_fun, points_matrix_q(:, 1), points_matrix_q(:, 2), 'UniformOutput', false);
          
        % The minimum distance for each point
        [distance_list, min_ind]    = cellfun(@min, distance_cell);
        
        % Resulting projections
        projected_points_matrix_q   = ellipse_samples(min_ind, :);
        
    %% Transforming back to original coordinate frame %%
        % They are no longer restricted to the first quarter
        projected_points_matrix_r   = sign_matrix .* projected_points_matrix_q;
        
        % Rotated and translated back to the original coordinate frame
        projected_points_matrix_t   = (ellipse_axes \ projected_points_matrix_r')';
        projected_points_matrix     = projected_points_matrix_t + ellipse_centre;
    
    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on    
            
            % Ellipse
            number_coord            = 1e3;
            ellipse_coord_matrix    = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, number_coord);    
            plot(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Ellipse');

            % Points
            scatter(points_matrix(:, 1), points_matrix(:, 2), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Points');

            % Projected points
            scatter(projected_points_matrix(:, 1), projected_points_matrix(:, 2), 'filled' ,'MarkerFaceColor', 'm', 'DisplayName', 'Projected points');

            
            % Axes
            axis equal
            xlabel('x [m]');
            ylabel('y [m]');

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            % Pause
            disp('The points have been projected. The script will finish and figure will close upon a key-press');
            pause();
            
            close(1);    
        end
    
end
    