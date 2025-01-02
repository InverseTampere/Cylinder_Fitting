% There can be 0 to 4 intersections between an ellipse and a circle, following from the solution to a quartic polynomial
% Theta is the counter-clockwise angle from the horizontal axis around the circle

function [intersection_matrix, number_intersections, theta_list] = Circle_Ellipse_Intersection(circle_centre, circle_radius, ellipse_centre, ellipse_radii, ellipse_axes, Print, Plot, Diagnostics)

    %% Coordinate system transformation %%
        % A circle-centered, ellipse-aligned coordinate system is used
        ellipse_centre_t = ellipse_centre - circle_centre;
        ellipse_centre_r = ellipse_centre_t / ellipse_axes;

    %% Intersections %%
        % The coefficients of the quartic
        s = 1/ellipse_radii(1)^2 - 1/ellipse_radii(2)^2;
        k = sum(ellipse_centre_r.^2 ./ ellipse_radii.^2) + circle_radius^2 / ellipse_radii(2)^2 - 1;

        a = s^2;
        b = -4*ellipse_centre_r(1)*s / ellipse_radii(1)^2;
        c = sum(4*ellipse_centre_r.^2 ./ ellipse_radii.^4) + 2*k*s;
        d = -4*ellipse_centre_r(1)*k / ellipse_radii(1)^2;
        e = k^2 - (2*ellipse_centre_r(2)*circle_radius / ellipse_radii(2)^2)^2;
    
        % Its unique roots
        margin = circle_radius * 1e-3;

        if s == 0           % If both radii are equal, s equals zero and the polynomial reduces to a quadratic one
            [x_intersection_r_list, ~, ~]   = Quadratic_Polynomial(c, d, e);
        else
            [x_intersection_r_list, ~, ~]   = Quartic_Polynomial(a, b, c, d, e, Diagnostics, Diagnostics, Diagnostics);
            [x_intersection_r_list, ~]      = Unique_Margin(x_intersection_r_list, margin);
        end
        
        % The y-coordinate of the intersections
        y_intersection_r_list = sqrt(circle_radius^2 - x_intersection_r_list.^2);
        
        % Complex intersections are removed
        imag_component          = abs(imag(y_intersection_r_list));
        x_intersection_r_list   = x_intersection_r_list(imag_component < margin);                       % A very small imaginary component may be present due to Matlab rounding
        y_intersection_r_list   = y_intersection_r_list(imag_component < margin);
                
        % The sign of the y-coordinate may be false, and is verified using the elliptical equation
        intersection_sign_matrix    = [x_intersection_r_list, y_intersection_r_list; x_intersection_r_list, -y_intersection_r_list];
        intersection_sign_matrix    = unique(intersection_sign_matrix, 'rows');                         % Removes duplicate rows in case of y being zero
        
        elliptical_equation_list    = sum((intersection_sign_matrix - ellipse_centre_r).^2 ./ ellipse_radii.^2, 2) - 1;
        sign_boolean_list           = abs(elliptical_equation_list) < margin;                           % The difference should be very close to 0
        
        intersection_matrix_r       = intersection_sign_matrix(sign_boolean_list, :);
        number_intersections        = size(intersection_matrix_r, 1);
        
    %% Transformation to original coordinate system %%
        % Resulting intersections in original coordinate frame
        intersection_matrix_t   = intersection_matrix_r * ellipse_axes;
        intersection_matrix     = intersection_matrix_t + circle_centre;
        
        % Angles around the circle
        theta_list                  = atan2(intersection_matrix(:, 2), intersection_matrix(:, 1));
        theta_list(theta_list < 0)  = 2*pi + theta_list(theta_list < 0);        % Ensures they are between 0 and 2pi

    %% Results %%
        %--% Printed messages %--%
        if Print == true
            fprintf('There are %g intersections: \n', number_intersections);
            
            for i = 1 : number_intersections
                fprintf('   [%.3g, %.3g] \n', intersection_matrix(i, 1), intersection_matrix(i, 2));
            end
        end
        
        %--% Figure %--%
        if Plot == true
            number_coord = 1e2;
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Circle
            [x_circle_list, y_circle_list] = Circle_Coordinates(circle_radius, circle_centre(1), circle_centre(2), number_coord);
            plot(x_circle_list, y_circle_list, 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Circle');

            % Ellipse
            ellipse_coord_matrix = Ellipse_Coordinate_Generator(ellipse_centre, ellipse_axes, ellipse_radii, number_coord);    
            plot(ellipse_coord_matrix(:, 1), ellipse_coord_matrix(:, 2), 'LineWidth', 2, 'color', 'r', 'DisplayName', 'Ellipse');

            % Intersections
            scatter(intersection_matrix(:, 1), intersection_matrix(:, 2), 'filled', 'MarkerFaceColor', 'm', 'DisplayName', 'Intersections');

            % Axes
            axis equal
            xlabel('x [m]');
            ylabel('y [m]');

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            % Pause message
            disp('The intersections have been found. The figure will close and script will end upon a key-press');
            pause();
            
            close(1);
        end
    
end
    