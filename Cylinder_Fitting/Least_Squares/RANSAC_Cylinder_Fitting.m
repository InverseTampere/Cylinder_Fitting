% Cylinder fitting through RANSAC of the least-squares cylinder fitting approach

function Cylinder_Geometry = RANSAC_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Scanner_Parameters, Fitting_Parameters, Statistical_Values, Plot)
    
    %% Structure inputs %%
        % Point cloud
        point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;
        number_points_list      = Point_Cloud_Coord.number_points_list;

        % Scanning parameters
        number_scanners         = Scanning_Parameters.number_scanners;
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;

        % Scanner parameters
        range_bias              = Scanner_Parameters.range_bias;

        % Statistical values
        Confidence_interval     = Statistical_Values.Confidence_interval;

    %% Manual inputs %%
        number_STD_step_size    = 0.2;          % [-] If a valid set is not found, the cylinder distance threshold is increased by this amount
        max_number_sets         = 1e2;          % [-] The number of sets may never exceed this amount
        Set_Size_Warning        = false;        % [true, false] A warning may be printed if the number of sets is limited
        
        LS_Diagnostics          = false;        % [true, false] Diagnostics of least-squares cylinder fitting
        Distr_Diagnostics       = false;        % [true, false] Generation of distributions based on the initial geometry estimate
        Bias_Diagnostics        = false;        % [true, false] Creates a plot of the expected bias in propagation error and cylinder distance

    %% Initial geometry and distribution estimate %%
        % Initial geometry fit comes from least-squares
        Init_Cylinder_Geometry  = Least_Squares_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, Confidence_interval, LS_Diagnostics, LS_Diagnostics);
        cylinder_centre_init    = Init_Cylinder_Geometry.Cylinder_centre;
        cylinder_radius_init    = Init_Cylinder_Geometry.Cylinder_radius;
        cylinder_direction_init = Init_Cylinder_Geometry.Cylinder_direction;

        % The range and radial uncertainty of each point
        [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(cylinder_centre_init, cylinder_direction_init, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
    
        % The 3D distributions
        alpha                       = 1 - Confidence_interval/100;          % The confidence interval is changed to alpha
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

    %% RANSAC %%
        % The number of required inliers
        number_points           = sum(number_points_list);
        min_number_points       = 5;                                    % As there are 5 infinite cylinder variables
        
        P_threshold             = Confidence_interval / 100;
        number_inlier_threshold = ceil(P_threshold * number_points);

        % Distances to the cylinder
        point_cloud_matrix  = vertcat(point_cloud_cell{:});
        [~, ~, omega_list]  = Point_to_Vector_Projection(point_cloud_matrix, cylinder_direction_init, cylinder_centre_init);
        delta_list          = abs(omega_list - cylinder_radius_init);
    
        [expected_delta_list, ~, ~, ~] = Cylinder_Expected_Bias(cylinder_radius_init, cylinder_centre_init, cylinder_direction_init, Point_Cloud_Distributions, Scanner_loc_cell, Bias_Diagnostics);
    
        % If a valid set is not found, the epsilon threshold is increased
        number_STD          = -number_STD_step_size;         % Note that this puts it at 0 for the first iteration
        max_number_inliers  = 0;

        while max_number_inliers <= number_inlier_threshold
            % The number of standard deviations is increased
            number_STD = number_STD + number_STD_step_size;
    
            % The error threshold follows from the expected values for delta
            expected_delta_STD  = std(expected_delta_list);
            mean_expected_delta = mean(expected_delta_list);
            delta_threshold     = mean_expected_delta + number_STD * expected_delta_STD;
            
            % Number of outliers    
            number_outliers = sum(delta_list > delta_threshold);
            P_outlier       = number_outliers / number_points;

            % If there are no outliers, all points are inliers by definition
            if P_outlier == 0
                max_number_inliers  = number_points;
                inlier_bool_cell    = {true(1, number_points)};
                max_ind             = 1;

                break            
            end
        
            % Number of random sets
            number_sets = ceil(log(1 - P_threshold) / log(1 - (1 - P_outlier)^min_number_points));

            % It may not exceed a maximum, which likely indicates a problem
            if number_sets > max_number_sets
                if Set_Size_Warning == true
                    warning('The number of sets %i exceeds the allowed maximum.', number_sets);
                end

                number_sets = max_number_sets;

            % Additionally if there are only outliers, the number of sets becomes -Inf and is therefore set at the maximum
            elseif P_outlier == 1
                number_sets = max_number_sets;
            end
        
            % Sets of point indices
            Set_Sampling_fun    = @(x) randsample(number_points, min_number_points);
            set_point_ind_cell  = arrayfun(Set_Sampling_fun, 1:number_sets, 'UniformOutput', false);
        
            % RANSAC fitting
            LS_Cylinder_Fitting_fun                 = @(set_point_ind) LS_Cylinder_Fitting(set_point_ind, point_cloud_matrix, Scanner_loc_cell, Fitting_Parameters, Confidence_interval, delta_threshold, LS_Diagnostics);
            [inlier_bool_cell, number_inliers_cell] = cellfun(LS_Cylinder_Fitting_fun, set_point_ind_cell, 'UniformOutput', false);
            
            % The set with the largest number of inliers
            number_inliers_list             = vertcat(number_inliers_cell{:});
            [max_number_inliers, max_ind]   = max(number_inliers_list);        
        end

        % The largest inlier set's point cloud
        inlier_bool = inlier_bool_cell{max_ind};    

        cumsum_points_list      = [0, cumsum(number_points_list)];
        inlier_point_cloud_cell = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's inliers are used as its point cloud
            first_ind   = cumsum_points_list(s) + 1;
            last_ind    = cumsum_points_list(s + 1);

            scanner_inlier_bool         = inlier_bool(first_ind : last_ind);
            inlier_point_cloud_cell{s}  = point_cloud_cell{s}(scanner_inlier_bool, :);
        end

        number_inlier_points_list   = cellfun(@length, inlier_point_cloud_cell);
        Inlier_Point_Cloud_Coord    = struct('point_cloud_cell', {inlier_point_cloud_cell}, 'number_points_list', number_inlier_points_list); 

        % Final least-squares fit over this inlier set
        Cylinder_Geometry = Least_Squares_Cylinder_Fitting(Inlier_Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, Confidence_interval, LS_Diagnostics, LS_Diagnostics);

    %% Plot %%
        if Plot == true
            % Printed message
            number_iterations = round(number_STD / number_STD_step_size + 1);       % Rounded to avoid rounding errors making it a float
            inlier_percentage = max_number_inliers / number_points * 100;
            fprintf('Fit achieved after %i iterations (%.2g sigma). Percentage of inliers is %.3g %% \n', number_iterations, number_STD, inlier_percentage);

            % The number of coordinates in the cylinder
            number_coord = 1e3;     

            % Figure
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid on
                
            % Initial cylinder surface
            cylinder_length_init = Init_Cylinder_Geometry.Cylinder_length;
            [init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length_init, cylinder_centre_init, cylinder_direction_init, number_coord);
            surf(init_cylinder_coord_x, init_cylinder_coord_y, init_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'c', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Initial cylinder');

            % Final cylinder surface
            [cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction] = deal(Cylinder_Geometry.Cylinder_radius, Cylinder_Geometry.Cylinder_length, Cylinder_Geometry.Cylinder_centre, Cylinder_Geometry.Cylinder_direction);
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Final cylinder');

            % The (inlier) point cloud
            scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Point cloud');
            scatter3(point_cloud_matrix(inlier_bool, 1), point_cloud_matrix(inlier_bool, 2), point_cloud_matrix(inlier_bool, 3), 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Inlier point cloud');

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

            disp('The RANSAC cylinder fitting script will finish and the plot will close when a key is pressed');
            pause();

            close(1);             
        end

    %% Local functions %%    
        % RANSAC fitting of a set with least-squares
        function [inlier_bool, number_inliers] = LS_Cylinder_Fitting(set_point_ind, point_cloud_matrix, Scanner_loc_cell, Fitting_Parameters, Confidence_interval, epsilon_threshold, LS_Diagnostics)
            % Input point cloud set - said to be from only one scanner due to the low number of points
            set_point_cloud         = point_cloud_matrix(set_point_ind, :);
            number_set_points       = length(set_point_ind);
            Set_Point_Cloud_Coord   = struct('point_cloud_cell', {{set_point_cloud}}, 'number_points_list', number_set_points);

            % Input scanning parameters only consider the first scanner
            Set_Scanning_Parameters = struct('Scanner_loc_cell', {Scanner_loc_cell(1)}, 'number_scanners', 1);

            % Filtering is always turned off as the set containts the minimum number of points
            Set_Fitting_Parameters              = Fitting_Parameters;
            Set_Fitting_Parameters.LS_Filtering = false;

            % Least-squares cylinder fitting
            LS_Cylinder_Geometry    = Least_Squares_Cylinder_Fitting(Set_Point_Cloud_Coord, Set_Scanning_Parameters, Set_Fitting_Parameters, Confidence_interval, LS_Diagnostics, LS_Diagnostics);
            cylinder_direction_LS   = LS_Cylinder_Geometry.Cylinder_direction;
            cylinder_centre_LS      = LS_Cylinder_Geometry.Cylinder_centre;
            cylinder_radius_LS      = LS_Cylinder_Geometry.Cylinder_radius;

            [~, ~, omega_list_LS]   = Point_to_Vector_Projection(point_cloud_matrix, cylinder_direction_LS, cylinder_centre_LS);
            epsilon_list_LS         = abs(omega_list_LS - cylinder_radius_LS);

            inlier_bool     = epsilon_list_LS < epsilon_threshold;
            number_inliers  = sum(inlier_bool);
        end
end
