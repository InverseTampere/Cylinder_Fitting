% The intersection between a plane and a sphere either does not exist, or is circular

function [circle_centre_matrix, circle_radius_list] = Plane_Sphere_Intersection(plane_point_matrix, plane_normal_vector_matrix, sphere_centre, sphere_radius)

    %% Manual inputs %%
        Plot    = false;     % [true, false] Shows the plane, sphere and the resulting circle (if it exists)
        
    %% Circle centre %%
        % Normalise the normal vectors
        plane_normal_vector_matrix = plane_normal_vector_matrix ./ sqrt(sum(plane_normal_vector_matrix.^2, 2));   

        % The distance from sphere centre to circle centre, along the plane normal vector
        alpha_list = dot(plane_point_matrix - sphere_centre, plane_normal_vector_matrix, 2);     
        
        % The resulting location of the circle centre
        circle_centre_matrix = alpha_list .* plane_normal_vector_matrix + sphere_centre;                                         

    %% Circle radius %%
        centre_vector_norms     = abs(alpha_list);                                      % The distances between the circle centres and the sphere centre
        circle_radius_list      = sqrt(sphere_radius^2 - centre_vector_norms.^2);       % The circle radius then follows from a simple right-angle triangle

    %% Valid intersection check %%
        % If there is no intersection, the radius is imaginary
        % It can however occur that there is a zero-radius intersection with a slight imaginary component
        imag_threshold          = 1e-2;                     % These are split using this threshold
        
        imag_components         = imag(circle_radius_list);
        invalid_intersections   = abs(imag_components) >  imag_threshold;      
        point_intersections     = abs(imag_components) > 0 & abs(imag_components) < imag_threshold;
        
        % NaN values are given to non-intersections
        circle_radius_list(invalid_intersections)       = NaN;
        circle_centre_matrix(invalid_intersections, :)  = NaN;
        
        % A radius of 0 is given to zero-radius intersections (i.e. the value without the imaginary component)
        circle_radius_list(point_intersections)         = 0;

    %% Plot %%
        if Plot == true
            % Number of points into which the spheres and circles are discretised
            n_discr     = 1e2;
            
            % Colour map for the planes
            num_planes  = size(plane_point_matrix, 1);
            
            plane_cmap  = cbrewer('qual', 'Set1', max(num_planes, 3));
            plane_cmap  = max(plane_cmap, 0);
            plane_cmap  = min(plane_cmap, 1);
            
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
            
            hold on
            grid on
    
            % The sphere
            [unit_sphere_x_matrix, unit_sphere_y_matrix, unit_sphere_z_matrix] = sphere(n_discr);
            surf(sphere_centre(1) + sphere_radius * unit_sphere_x_matrix, sphere_centre(2) + sphere_radius * unit_sphere_y_matrix, sphere_centre(3) + sphere_radius * unit_sphere_z_matrix, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'DisplayName', 'Sphere');
            
            % The results for each given plane
            for p = 1 : num_planes
                plane_colour = plane_cmap(p, :);
                
                % The normal vector to the plane
                plane_point         = plane_point_matrix(p, :);
                plane_normal_vector = plane_normal_vector_matrix(p, :);
                plane_normal_vector = -plane_normal_vector * sphere_radius / 2;  % Scaled and reversed
                
                pl_normal = plot3(plane_point(1) + [0, plane_normal_vector(1)], plane_point(2) + [0, plane_normal_vector(2)], plane_point(3) + [0, plane_normal_vector(3)], 'color', plane_colour, 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', 'Plane normal');
            
                % The circle
                circle_centre = circle_centre_matrix(p, :);
                circle_radius = circle_radius_list(p);
                
                [a_circle_list, b_circle_list]  = Circle_Coordinates(circle_radius, 0, 0, n_discr);                     % 2D coordinates
                circle_matrix_r                 = [a_circle_list', b_circle_list', zeros(n_discr, 1)];                  % The 'height' is zero in the circle's own coordinate frame

                [~, vector_basis, ~]    = Vector_Based_Rotation(circle_centre, plane_normal_vector, circle_centre);
                circle_matrix_t         = (vector_basis' * circle_matrix_r')';                                          % Rotation by inverse vector basis
                circle_matrix           = circle_matrix_t + circle_centre;                                              % Translated

                pl_circle = plot3(circle_matrix(:, 1), circle_matrix(:, 2), circle_matrix(:, 3), 'LineWidth', 1, 'color', plane_colour, 'DisplayName', 'Intersection circle');
                
                % The plane
                sphere_data_matrix  = repmat(sphere_radius * [-1; 1], [1, 3]);
                plane_corner_matrix = Plane_Corner_Points(plane_normal_vector, plane_point, [plane_point; sphere_data_matrix]);
        
                pa_plane = patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), plane_colour, 'FaceAlpha', 0.5, 'DisplayName', 'Plane');

                if p > 1
                    pl_normal.HandleVisibility  = 'Off';
                    pl_circle.HandleVisibility  = 'Off';
                    pa_plane.HandleVisibility   = 'Off';
                end
            end
            
            % Axes
            xlabel('x');
            ylabel('y');
            zlabel('z');
            
            axis equal
            
            view(45, 45);
            
            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    

            % Time to look at the figure
            disp('The figure is closed when a button is pressed, and the script continues.');
            pause();
            close(1);
        end
        
end