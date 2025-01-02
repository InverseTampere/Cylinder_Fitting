% The infinite cylinder (length is not parametrised) is fitted through least-squares by minimising the error function arising from least-squares,
% as a function of the azimuth/elevation angle of the plane normal to the axis

% Note that the circle centre is relative to the point cloud centroid and normal to the cylinder axis
% If weighting is not desired, an empty input may be given

function [circle_centre, cylinder_centre, cylinder_radius, cylinder_direction] = Least_Squares_Inf_Cylinder_Fitting(point_matrix, weight_list, Plot)
    
    %% Manual inputs %%
        convergence_threshold   = 1e-2;         % Convergence threshold for power iteration, as first estimate for the optimiser
        
    %% The centered point cloud %%
        point_centroid  = mean(point_matrix, 1);        
        point_matrix_c  = point_matrix - point_centroid;
        num_points      = size(point_matrix_c, 1);
        
    %% Initial vector angles estimate %%
        % Initial estimate of the spherical angles from the principal component
        covariance_matrix = cov(point_matrix_c);
        [init_cyl_dir, ~] = Power_Iteration(covariance_matrix, convergence_threshold);
    
        init_azim = atan(init_cyl_dir(2) / init_cyl_dir(1));
        init_elev = atan(init_cyl_dir(3) / sqrt(init_cyl_dir(1)^2 + init_cyl_dir(2)^2));
        
        init_spherical_angles = [init_azim, init_elev];

        % Ensure they are within the bounds
        init_spherical_angles = sign(init_spherical_angles) .* init_spherical_angles;                                   % Ensures the angles are positive
        init_spherical_angles(init_spherical_angles > pi) = init_spherical_angles(init_spherical_angles > pi) - pi;     % Ensures they are between [0, pi]

    %% Minimisation of the least squares error %%
        % The weights
        if isempty(weight_list)     
            weight_list = ones(num_points, 1);                  % If weights are not given, they are uniform
        else
            weight_list = weight_list / mean(weight_list);      % To ensure that the sum is a reasonable value
        end

        % The bounds of the spherical angles
        lower_bounds    = zeros(1, 2);
        upper_bounds    = pi*ones(1, 2);
        
        % The optimisation options
        Options = optimoptions('fmincon', 'optimalitytolerance', 1e-6, 'steptolerance', 1e-6, 'Display', 'Off');
                  
        % Optimisation
        Least_Squares_fun   = @(spherical_angles) Least_Squares_Error(spherical_angles);
        spherical_angles    = fmincon(Least_Squares_fun, init_spherical_angles, [], [], [], [], lower_bounds, upper_bounds, [], Options);
        
    %% The resulting geometry %%
        % Circle geometry and cylinder direction
        [circle_centre_c, cylinder_radius, cylinder_direction, vector_basis, ~] = Least_Squares_Circle_Fitting(spherical_angles);

        % The circle centre is moved back to account for the point cloud centering
        rotated_point_centroid  = (vector_basis * point_centroid')';
        
        proj_point_centroid     = rotated_point_centroid(1 : 2);
        height_point_centroid   = rotated_point_centroid(3);
    
        circle_centre   = circle_centre_c + proj_point_centroid;
        
        % Cylinder centre     
        cylinder_centre = (vector_basis \ [circle_centre, height_point_centroid]')';
        
    %% Plot %%
        % Plot showing the point cloud and fitted cylinder
        if Plot == true
            % The number of coordinates in the cylinder
            number_coord = 1e3;     
            
            % Estimate of the length
            [~, delta_list]     = Point_to_Vector_Projection(point_matrix_c, cylinder_direction, cylinder_centre);
            cylinder_length     = 2*max(delta_list);            
            
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % The initial cylinder direction
            plot3(point_centroid(1) + cylinder_length/2 * [0, init_cyl_dir(1)], point_centroid(2) + cylinder_length/2 * [0, init_cyl_dir(2)], point_centroid(3) + cylinder_length/2 * [0, init_cyl_dir(3)], 'LineWidth', 2, 'color', 'b', 'DisplayName', 'Initial direction');
            
            % The fitted cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'g', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Least-squares fitted infinite cylinder');

            % The point cloud
            scatter3(point_matrix(:, 1), point_matrix(:, 2), point_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Point cloud');
            
            % Aspect ratio
            axis equal

            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')
            
            % Viewing angle
            view(45, 45)

            % Legend
            legend('show', 'location', 'northoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            disp('The script will continue and the plot will close when a key is pressed');
            pause();

            close(1);             
        end
        
    %% Local functions %%
        % Least-squares circle fitting
        function [circle_centre, circle_radius, cylinder_direction, vector_basis, proj_point_matrix] = Least_Squares_Circle_Fitting(spherical_angles)
            % The spherical angles
            azim = spherical_angles(1);
            elev = spherical_angles(2);
            
            % The resulting cylinder and plane vectors
            cylinder_direction = [cos(elev)*cos(azim), cos(elev)*sin(azim), sin(elev)];
            
            num_dim = length(cylinder_direction);
            origin  = zeros(1, num_dim);
            [~, vector_basis, plane_vectors] = Vector_Based_Rotation(origin, cylinder_direction, origin);

            % The projected point cloud
            proj_point_matrix = (plane_vectors * point_matrix_c')';
            
            % Least-squares circle fitting
            A = [proj_point_matrix, ones(num_points, 1)];
            B = -sum(proj_point_matrix.^2, 2);
            
            circle_coefficients = A \ B;
            circle_centre       = -1/2 * circle_coefficients(1:2)'; 
            circle_radius       = sqrt(sum(circle_centre.^2) - circle_coefficients(3));
        end
    
        % The least-squares error
        function avg_geometry_error = Least_Squares_Error(spherical_angles)
            % The least-squares geometry
            [LS_circle_centre, LS_circle_radius, ~, ~, proj_point_matrix] = Least_Squares_Circle_Fitting(spherical_angles);
            
            % The geometry error
            centre_dist_list    = sqrt(sum((proj_point_matrix - LS_circle_centre).^2, 2));
            abs_error_list      = abs(centre_dist_list - LS_circle_radius);
            abs_error_list_w    = weight_list .* abs_error_list;                    % Weighted
            avg_geometry_error  = sum(abs_error_list_w) / num_points;               % Normalised by the number of points
        end

end