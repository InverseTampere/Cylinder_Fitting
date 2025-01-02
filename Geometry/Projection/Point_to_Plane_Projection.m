% This script projects a 3D coordinate matrix onto a 2D plane that lies on the given point.
% The plane is normal to the 3rd axis of the given vector basis
% The projected coordinates are in the original coordinate frame, the plane coordinates are according to the given vector basis and plane point

function [proj_coord_matrix, plane_coord_matrix, height_list] = Point_to_Plane_Projection(coord_matrix, vector_basis, plane_point, Plot)
    
    %% Coordinate matrix projection %%
        % The coordinate matrix is translated and rotated to the vector basis
        transl_coord_matrix = coord_matrix - plane_point;
        rot_coord_matrix    = (vector_basis * transl_coord_matrix')';
        
        % The third component (in the vector direction) is separated
        [num_points, num_dim]   = size(coord_matrix);
        plane_coord_matrix      = [rot_coord_matrix(:, 1:2), zeros(num_points, 1)];
        height_list             = rot_coord_matrix(:, 3);
        
        % The coordinates in the original coordinate frame
        proj_coord_matrix_t     = (vector_basis \ plane_coord_matrix')';
        proj_coord_matrix       = proj_coord_matrix_t + plane_point;
        
    %% Plot %%
        if Plot == true
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])

            %--% The original coordinate frame %--%
            subplot(1, 2, 1)
            hold on
            grid on
            
            % The original coordinate matrix
            scatter3(coord_matrix(:, 1), coord_matrix(:, 2), coord_matrix(:, 3), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Original coordinates');            

            % The projected coordinate matrix
            scatter3(proj_coord_matrix(:, 1), proj_coord_matrix(:, 2), proj_coord_matrix(:, 3), 'filled', 'MarkerFaceColor', 'g', 'DisplayName', 'Projected coordinates');            
            
            % The scaled plane normal vector
            data_ampl           = max(coord_matrix) - min(coord_matrix);
            lambda              = sqrt(sum(data_ampl.^2));
            
            plane_normal_vector = lambda * vector_basis(num_dim, :);
            
            plot3(plane_point(1) + [0, plane_normal_vector(1)], plane_point(2) + [0, plane_normal_vector(2)], plane_point(3) + [0, plane_normal_vector(3)], 'LineWidth', 2, 'color', 'm', 'DisplayName', 'Plane normal vector');
            
            % Plane
            plane_corner_matrix = Plane_Corner_Points(plane_normal_vector, plane_point, [proj_coord_matrix; coord_matrix]);
            patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'm', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');               
            
            % The axes
            xlabel('x');
            ylabel('y');
            zlabel('z');

            % The aspect ratio is made equal
            axis equal

            % Legend
            legend('show');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % The angle from which the plot is viewed
            view(45, 45);
            
            %--% The rotated coordinate frame %--%
            subplot(1, 2, 2)
            hold on
            grid on
            
            % The rotated coordinate matrix
            scatter3(rot_coord_matrix(:, 1), rot_coord_matrix(:, 2), rot_coord_matrix(:, 3), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Original coordinates');            

            % The projected coordinate matrix
            scatter3(plane_coord_matrix(:, 1), plane_coord_matrix(:, 2), plane_coord_matrix(:, 3), 'filled', 'MarkerFaceColor', 'g', 'DisplayName', 'Projected coordinates');            
            
            % The scaled plane normal vector
            plot3([0, 0], [0, 0], [0, lambda], 'LineWidth', 2, 'color', 'm', 'DisplayName', 'Plane normal vector');
            
            % Plane
            z_axis  = [0, 0, 1];
            origin  = zeros(1, num_dim);
            plane_corner_matrix = Plane_Corner_Points(z_axis, origin, [plane_coord_matrix; rot_coord_matrix]);
            patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), 'm', 'FaceAlpha', 0.25, 'DisplayName', 'Plane');               
            
            % The axes
            xlabel('x_v');
            ylabel('y_v');
            zlabel('z_v');

            % The aspect ratio is made equal
            axis equal

            % Legend
            legend('show');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % The angle from which the plot is viewed
            view(45, 45);

            % Pause
            disp('The plot will close and script continue if a key is pressed');
            pause();
            
            close(1)
        end
    
end