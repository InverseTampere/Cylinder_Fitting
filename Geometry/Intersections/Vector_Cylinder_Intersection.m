% This script computes the intersection points between line vectors and a cylinder
% In the given matrices, the columns stand for the [x, y, z] coordinates respectively
% Note that NaN coordinates are returned if an intersection does not exist

% The geometry is viewed as a cylinder with flat end caps. Intersections can occur on these end caps as well, 
% however are marked as true/false such that they can be easily removed if these aren't desired

function [intersection_point_matrix, cylinder_intersects, top_end_cap_intersects, bot_end_cap_intersects, incidence_angle_list] = Vector_Cylinder_Intersection(line_origin_matrix, line_vec_matrix, cyl_radius, cyl_centre, cyl_length, cyl_direction)
    
    %% Inputs %%
        Diagnostics = false;        % [true, false] Shows the intersections, their surface normals and the resulting incidence angles

    %% The normalised direction vectors %% 
        % The direction of the cylinder
        cyl_direction_n     = cyl_direction / norm(cyl_direction);

        % The normalised line vectors
        line_vec_norms      = sqrt(sum(line_vec_matrix.^2, 2));
        line_vec_matrix_n   = line_vec_matrix ./ line_vec_norms;
                
    %% The intersection point of each vector %%
        number_vectors = size(line_origin_matrix, 1);
        
        % Intersections can occur on the cylinder or top/bottom endcaps (as seen from the cylinder)
        cylinder_intersects         = false(number_vectors, 1);
        top_end_cap_intersects      = false(number_vectors, 1);
        bot_end_cap_intersects      = false(number_vectors, 1);
        
        % The coordinates of the intersection points are NaN if one isn't found
        intersection_point_matrix   = NaN(number_vectors, 3);
        
        % The surface normal vectors
        surface_normal_matrix       = NaN(number_vectors, 3);
        
        % The incidence angles
        incidence_angle_list        = NaN(number_vectors, 1);
        
        for v = 1 : number_vectors
            % This vector's start point and direction
            line_origin = line_origin_matrix(v, :);
            line_vec_n  = line_vec_matrix_n(v, :);
                        
            % The origin of the cylinder is translated w.r.t. the origin of this line
            cyl_centre_t    = cyl_centre - line_origin;
            
            %--% End cap intersections %--%
            top_cap     = cyl_centre_t + cyl_direction_n * cyl_length / 2;
            bot_cap     = cyl_centre_t - cyl_direction_n * cyl_length / 2;
            
            lambda_top  = dot(cyl_direction_n, top_cap) / dot(cyl_direction_n, line_vec_n);
            lambda_bot  = dot(cyl_direction_n, bot_cap) / dot(cyl_direction_n, line_vec_n);
            
            % The difference between the intersects and end caps
            difference_top  = line_vec_n * lambda_top - top_cap;
            difference_bot  = line_vec_n * lambda_bot - bot_cap;
            
            % Intersections only exist if the end cap planes are intersected, and then if they are intersected within the radius of the caps
            planar_intersection_top = abs(dot(difference_top, cyl_direction_n)) < 1e-6;         % Should be zero, but non-zero value to account for Matlab rounding
            planar_intersection_bot = abs(dot(difference_bot, cyl_direction_n)) < 1e-6;
            
            if ~planar_intersection_top || norm(difference_top) > cyl_radius                    % Additionally, an intersect also has to be within the radius
                lambda_top = NaN;                                                               % If it doesn't intersect, a NaN value is given
            end

            if ~planar_intersection_bot || norm(difference_bot) > cyl_radius
                lambda_bot = NaN;
            end
            
            %--% Cylinder intersections %--%
            % The cross product between the two direction vectors
            vec_cross   = cross(line_vec_n, cyl_direction_n);
            
            % The determinant is computed to check for the number of intersections
            determinant = cyl_radius^2 * dot(vec_cross, vec_cross) - dot(cyl_centre_t, vec_cross)^2;
            
            if determinant > 0         % Only in this case a solution exists
                % Parts of the quadratic equation
                numerator       = dot(vec_cross, cross(cyl_centre_t, cyl_direction_n));
                denominator     = dot(vec_cross, vec_cross);

                lambda_list_cyl = (numerator + sqrt(determinant) * [-1, 1]) / denominator;
                
                % Check if the intersection point lies within the finite cylinder length
                cyl_direction_matrix_n      = cyl_direction_n' .* ones(1, 2);       % Duplicates for the two different distances
                cyl_axis_lambda_matrix_n    = line_vec_n' .* lambda_list_cyl - cyl_centre_t';
                cylinder_axis_lambda_list   = dot(cyl_direction_matrix_n, cyl_axis_lambda_matrix_n);

                cylinder_axis_distance_list = abs(cylinder_axis_lambda_list);
                
                lambda_list_cyl(cylinder_axis_distance_list > cyl_length / 2) = NaN;    % Note that half the length is used as it is assumed to be equal away from the centre
            
            else
                lambda_list_cyl = [];
            end
            
            %--% The closest intersect %--%
            % The list of end cap and cylinder lambdas
            lambda_list = [lambda_top, lambda_bot, lambda_list_cyl];
            
            % The first intersection point is taken (as the other one lies on the opposite side of the cylinder)
            distance_list       = abs(lambda_list);
            closest_intersect   = find(distance_list == min(distance_list), 1);

            if isempty(closest_intersect)       % No intersection point meets the requirements above
                lambda = NaN;
            else
                lambda = lambda_list(closest_intersect);
                
                % If the closest intersect is 1 or 2, the end caps were intersected
                if      closest_intersect == 1
                    top_end_cap_intersects(v) = true;
                elseif  closest_intersect == 2
                    bot_end_cap_intersects(v) = true;
                elseif  closest_intersect >= 3
                    cylinder_intersects(v) = true;
                end
            end

            % The resulting intersection point
            intersection_point = lambda * line_vec_n + line_origin;       % Note the translation       

            intersection_point_matrix(v, :) = intersection_point;
            
            %--% The incidence angle %--%
            % The vector normal to the surface at the intersection point
            if      closest_intersect == 1      % The top cap intersected
                surface_normal = cyl_direction_n;
                
            elseif  closest_intersect == 2      % The bot cap intersected
                surface_normal = -cyl_direction_n;
                
            elseif  closest_intersect >= 3      % The cylinder intersected
                % The intersection point is projected onto the cylinder axis
                [projected_intersection_point, ~] = Point_to_Vector_Projection(intersection_point, cyl_direction_n, cyl_centre);
                
                % The surface normal is the vector from the projected intersection to the actual intersection
                surface_normal  = intersection_point - projected_intersection_point;
                surface_normal  = surface_normal / norm(surface_normal); 
            else
                surface_normal  = NaN(1, 3);
            end
            
            surface_normal_matrix(v, :) = surface_normal;
            
            % The incidence angle
            incidence_angle         = acos(abs(dot(surface_normal, line_vec_n)));
            incidence_angle         = min(incidence_angle, pi/2 - incidence_angle);
            incidence_angle_list(v) = incidence_angle;
        end

    %% Diagnostics plot %%
    if Diagnostics == true
        % Incidence angle colourmap
        num_colours = 1e3;
        alpha_cmap  = cbrewer('seq', 'YlOrRd', num_colours);
        alpha_cmap  = max(alpha_cmap, 0);
        alpha_cmap  = min(alpha_cmap, 1);
        
        total_intersects    = logical(cylinder_intersects + top_end_cap_intersects + bot_end_cap_intersects);
        alpha_ind           = round((num_colours - 1) * incidence_angle_list(total_intersects)/(pi/2)) + 1;        % Note that the maximum incidence angle is pi/2

        % The intersection points
        intersection_point_matrix_total = intersection_point_matrix(total_intersects, :);
        
        intersection_point_matrix_cyl   = intersection_point_matrix(cylinder_intersects, :);
        intersection_point_matrix_top   = intersection_point_matrix(top_end_cap_intersects, :);
        intersection_point_matrix_bot   = intersection_point_matrix(bot_end_cap_intersects, :);
        
        %--% Plot %--%
        figure(1)
        % Size and white background
        set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
        set(gcf, 'color', [1, 1, 1])    

        for s = 1 : 2
            subplot(1, 2, s)
            hold on
            grid on

            if s == 1
                % The cylinder
                number_coord = 1e3;
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cyl_radius, cyl_length, cyl_centre, cyl_direction_n, number_coord);
                surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', 'Cylinder');

                intersections = find(total_intersects == true)';
                for v = intersections
                    % The intersecting laser beams
                    intersection_point  = intersection_point_matrix(v, :);
                    line_origin         = line_origin_matrix(v, :);

                    pl_int = plot3([line_origin(1), intersection_point(1)], [line_origin(2), intersection_point(2)], [line_origin(3), intersection_point(3)], 'color', 'r', 'LineWidth', 1, 'DisplayName', 'Laser beams');

                    % Their surface normals
                    surface_normal  = surface_normal_matrix(v, :);
                    pl_surf_norm    = plot3(intersection_point(1) + cyl_radius*[0, surface_normal(1)], intersection_point(2) + cyl_radius*[0, surface_normal(2)], intersection_point(3) + cyl_radius*[0, surface_normal(3)], 'color', 'b', 'LineWidth', 1, 'DisplayName', 'Surface normal');
                    
                    if v > 1
                        pl_int.HandleVisibility = 'Off';
                        pl_surf_norm.HandleVisibility = 'Off';
                    end
                end

                % The intersection points, constant colour
                scatter3(intersection_point_matrix_cyl(:, 1), intersection_point_matrix_cyl(:, 2), intersection_point_matrix_cyl(:, 3), 100, 'c', 'Marker', '.', 'DisplayName', 'Cylinder intersections');
                scatter3(intersection_point_matrix_top(:, 1), intersection_point_matrix_top(:, 2), intersection_point_matrix_top(:, 3), 100, 'm', 'Marker', '.', 'DisplayName', 'Top cap intersections');
                scatter3(intersection_point_matrix_bot(:, 1), intersection_point_matrix_bot(:, 2), intersection_point_matrix_bot(:, 3), 100, 'b', 'Marker', '.', 'DisplayName', 'Bottom cap intersections');

            else
                % The intersection points
                scatter3(intersection_point_matrix_total(:, 1), intersection_point_matrix_total(:, 2), intersection_point_matrix_total(:, 3), 100, alpha_cmap(alpha_ind, :), 'Marker', '.', 'DisplayName', 'All intersections');

                % Colorbar
                colormap(alpha_cmap);

                cb = colorbar;
                shading interp
                clim([0, pi/2])
                ylabel(cb, '\alpha [rad]');
                cb.FontSize = 25;        
            end

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')

            x_lim = [min(cylinder_coord_x(:)) - cyl_radius, max(cylinder_coord_x(:)) + cyl_radius];
            y_lim = [min(cylinder_coord_y(:)) - cyl_radius, max(cylinder_coord_y(:)) + cyl_radius];
            z_lim = [min(cylinder_coord_z(:)) - cyl_radius, max(cylinder_coord_z(:)) + cyl_radius];

            xlim(x_lim);
            ylim(y_lim);
            zlim(z_lim);

            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
        end

        disp('The script will continue and the plot will close when a key is pressed')
        pause()

        close(1)
    end
end