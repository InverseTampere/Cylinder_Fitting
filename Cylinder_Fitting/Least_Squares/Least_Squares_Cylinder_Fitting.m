% The cylinder is fitted through the least-squares procedure implemented in TreeQSM
% Note that the circle centre is relative to the point cloud centroid, and normal to the cylinder axis

function Cylinder_Geometry = Least_Squares_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, confidence_interval, Plot, Diagnostics)
    
    %% Structure inputs %%
        % Point cloud
        point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;
        number_points_list      = Point_Cloud_Coord.number_points_list;

        % Scanning parameters
        number_scanners         = Scanning_Parameters.number_scanners;
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;

        % Fitting parameters
        LS_Filtering            = Fitting_Parameters.LS_Filtering;
        point_weight_list       = Fitting_Parameters.point_weights_list;

    %% Manual inputs %%
        alpha_max_flip          = 20;           % [deg] The location of points with an incidence angle below this value is compared to that of the first intersections
                                                %       to see whether the cylinder is in front of the point cloud (and thus needs to be flipped) or not

    %% Outlier filtering %%
        % When uncertainty is of similar or greater magnitude than the cylinder length, the principal axis of the point cloud becomes tilted away from the true axis
        % To circumvent this, outliers are removed that are outside of the given confidence interval
        if LS_Filtering == true
            % Filtering bounds are symmetric
            percentile_LB   = (100 - confidence_interval)/2;
            percentile_UB   = 100 - percentile_LB;
            
            % Points are filtered separately for each scanner to prevent bias
            outlier_boolean_cell = cell(1, number_scanners);

            for s = 1 : number_scanners
                % The points are rotated s.t. the vector from scanner to centroid is on the z-axis
                scanner_point_cloud_matrix      = point_cloud_cell{s};
                number_points                   = number_points_list(s);

                scanner_point_cloud_centroid    = mean(scanner_point_cloud_matrix, 1);
                scanner_loc                     = Scanner_loc_cell{s};
                scanner_vector                  = scanner_point_cloud_centroid - scanner_loc;

                [scanner_point_cloud_matrix_r, ~, ~] = Vector_Based_Rotation(scanner_point_cloud_matrix, scanner_vector, scanner_point_cloud_centroid);

                % Outliers in each dimension
                num_dim         = length(scanner_loc);
                dim_outlier_boolean = false(number_points, num_dim);

                for d = 1 : num_dim
                    point_cloud_dim_list        = scanner_point_cloud_matrix_r(:, d);
                    dim_outlier_boolean         = isoutlier(point_cloud_dim_list, 'Percentiles', [percentile_LB, percentile_UB]);
                    dim_outlier_boolean(:, d)   = dim_outlier_boolean;
                end

                % Points are removed that are outliers in any dimension
                outlier_boolean         = logical(sum(dim_outlier_boolean, 2));
                outlier_boolean_cell{s} = outlier_boolean;
            end
        else
            number_points   = sum(number_points_list);
            outlier_boolean = false(number_points, 1);
        end

        % The points for each scanner are aggregated into one matrix and then outliers are filtered out
        point_cloud_matrix      = vertcat(point_cloud_cell{:});

        point_cloud_matrix_f    = point_cloud_matrix(~outlier_boolean, :);
        point_weight_list_f     = point_weight_list(~outlier_boolean);

    %% Cylinder fitting %%
        %--% Initial estimate %--%
        % The input for least-squares cylinder fitting comes from an iterative least-squares circle fitting procedure
        [~, cylinder_centre_init, cylinder_radius_init, cylinder_direction_init] = Least_Squares_Inf_Cylinder_Fitting(point_cloud_matrix_f, point_weight_list_f, Diagnostics);

        % If there is only one scanner, the cylinder may be placed in front of the point cloud
        if number_scanners == 1
            % An estimate of the cylinder length
            [~, delta_list, ~]      = Point_to_Vector_Projection(point_cloud_matrix_f, cylinder_direction_init, cylinder_centre_init);
            cylinder_length_init    = max(delta_list) - min(delta_list);

            % Intersections of the laser beams to the cylinder
            number_points       = size(point_cloud_matrix_f, 1);
            scanner_loc         = Scanner_loc_cell{1};
            scanner_loc_rep     = repmat(scanner_loc, number_points, 1);                % The scanner location is the starting point of all vectors

            point_vector_matrix = point_cloud_matrix_f - scanner_loc;                   % Vectors from the scanner to each point
            point_norm_list     = sqrt(sum(point_vector_matrix.^2, 2));
            point_vector_matrix = point_vector_matrix ./ point_norm_list;

            [intersection_point_matrix, ~, ~, ~, incidence_angle_list] = Vector_Cylinder_Intersection(scanner_loc_rep, point_vector_matrix, cylinder_radius_init, cylinder_centre_init, cylinder_length_init, cylinder_direction_init);

            % Distances from the points to the intersects
            intersect_vector_matrix = intersection_point_matrix - scanner_loc;
            intersect_norm_list     = sqrt(sum(intersect_vector_matrix.^2, 2));
            epsilon_list            = point_norm_list - intersect_norm_list;

            % At low incidence angles, the distance from the initial intersect to the point is generally greater than the radius if the location is flipped
            alpha_max_flip  = deg2rad(alpha_max_flip);
            flip_bool       = incidence_angle_list < alpha_max_flip;
            
            epsilon_flip_list = epsilon_list(flip_bool);
            
            if mean(epsilon_flip_list) > cylinder_radius_init
                % The cylinder centre is displaced by the diameter to place it behind the point cloud
                cylinder_centre_vector  = (cylinder_centre_init - scanner_loc) / norm(cylinder_centre_init - scanner_loc);
                cylinder_centre_init    = cylinder_centre_init + cylinder_centre_vector * 2*cylinder_radius_init;
            end
        end
        
        % A structure is created for the geometry
        Initial_Inf_Cylinder_Geometry = struct('radius', cylinder_radius_init, 'axis', cylinder_direction_init, 'start', cylinder_centre_init);

        %--% Final estimate %--%
        % The geometry estimate is computed using the least-squares cylinder fitting algorithm implemented in TreeQSM
        TreeQSM_Cylinder_Geometry                                               = least_squares_cylinder_2(point_cloud_matrix_f, Initial_Inf_Cylinder_Geometry, point_weight_list_f, []);
        [cylinder_radius, cylinder_direction, cylinder_bottom, cylinder_length] = deal(TreeQSM_Cylinder_Geometry.radius, TreeQSM_Cylinder_Geometry.axis, TreeQSM_Cylinder_Geometry.start, TreeQSM_Cylinder_Geometry.length);        
                
        % The least-squares algorithm computes the location of the bottom of the cylinder, which is translated to the centre used in this script
        cylinder_centre = cylinder_bottom + cylinder_length/2 * cylinder_direction;
        
        % The top location follows from the length, vector and bottom location
        cylinder_top = cylinder_bottom + cylinder_length * cylinder_direction;

        % Sometimes the TreeQSM algorithm produces significant errors
        % Reasonably, the radius and length cannot be greater than twice the largest distance between the centroid and any given point
        % Additionally, the radius cannot be negative
        point_cloud_centroid    = mean(point_cloud_matrix_f, 1);
        norm_list               = sqrt(sum((point_cloud_matrix_f - point_cloud_centroid).^2, 2));

        if cylinder_radius > 2 * max(norm_list) || cylinder_length > 2 * max(norm_list) || cylinder_radius <= 0
            % The initial values are used for the cross-section geometry
            cylinder_centre     = cylinder_centre_init;
            cylinder_radius     = cylinder_radius_init;
            cylinder_direction  = cylinder_direction_init;

            % The length and top and bottom locations are estimated using the maximum norm
            cylinder_length     = 2*max(norm_list);
            cylinder_bottom     = cylinder_centre - cylinder_direction * cylinder_length/2;
            cylinder_top        = cylinder_centre + cylinder_direction * cylinder_length/2;
        end

        % The final geometry is saved in a structure
        Cylinder_Geometry = struct('Cylinder_centre', cylinder_centre, 'Cylinder_direction', cylinder_direction, 'Cylinder_radius', cylinder_radius, 'Cylinder_length', cylinder_length, 'Top_loc', cylinder_top, 'Bottom_loc', cylinder_bottom);

    %% Plot %%    
        % Plot showing the point cloud and fitted cylinder        
        if Plot == true
            % The number of coordinates in the cylinder
            number_coord = 1e3;     

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % Initial (infinite) cylinder surface
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length, cylinder_centre_init, cylinder_direction_init, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'c', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial infinite cylinder');

            % Final cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Final cylinder');

            % The point cloud
            point_cloud_matrix = vertcat(point_cloud_cell{:});
            scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Point cloud');
            scatter3(point_cloud_matrix_f(:, 1), point_cloud_matrix_f(:, 2), point_cloud_matrix_f(:, 3), 'filled', 'MarkerFaceColor', 'g', 'DisplayName', 'Filtered point cloud');

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

            disp('The least-squares cylinder fitting script will finish and the plot will close when a key is pressed');
            pause();

            close(1);             
        end
end