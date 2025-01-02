% This script computes the expected distance from the point cloud to the circle, given their uncertainty and incidence angles
% The power n of the distance can be taken such that it is more sensitive to outliers, and it is equivalent to the nth moment of distance

function [Average_expected_distance, Expected_distance_cell] = Expected_Circle_Distance(distance_moment, Scanner_loc_cell, x_points_cell, y_points_cell, sigma_radial_cell, sigma_range_cell, range_bias, laser_radius_cell, Tree_centre_x, Tree_centre_y, Tree_radius, Point_Mirroring, Point_Duplication, Plot, Print)

    %% Inputs %%
        number_steps    = 1e2;  % Number of steps into which the integral is discretised as a Riemann sum in either direction
                                % 1e2 is recommended
        number_STD      = 3;    % The number of standard deviations in the propagation direction that are considered
                                
    %% Point cloud preprocessing %%
         %--% % Range bias compensation %--%
        [x_points_cell, y_points_cell] = Range_Bias_Compensation(x_points_cell, y_points_cell, range_bias, Scanner_loc_cell);
                       
        %--% Point mirroring %--%
        number_scanners = length(Scanner_loc_cell);

        if Point_Mirroring == true    
            x_int_points_cell   = cell(1, number_scanners);
            y_int_points_cell   = cell(1, number_scanners);

            for s = 1 : number_scanners
                % This scanner's location and point cloud
                Scanner_loc_x   = Scanner_loc_cell{s}(1);
                Scanner_loc_y   = Scanner_loc_cell{s}(2);

                x_points_list   = x_points_cell{s};
                y_points_list   = y_points_cell{s};

                number_beams    = length(x_points_list);

                %--% Intersection points %--%
                % For beams that intersect the circle, the intersection point is used  
                x_start_list        = Scanner_loc_x * ones(number_beams, 1);
                y_start_list        = Scanner_loc_y * ones(number_beams, 1);

                vector_beams_x_list = x_points_list - x_start_list;
                vector_beams_y_list = y_points_list - y_start_list;

                [x_i_list, y_i_list, ~, intersection_list]  = Vector_Circle_Intersection(x_start_list, y_start_list, vector_beams_x_list, vector_beams_y_list, Tree_centre_x, Tree_centre_y, Tree_radius);

                % If no intersection point is used, the closest point on the circle to the point is used
                missed_beams = intersection_list == false;

                % The intersection points are the closest point on the circle to the point cloud, i.e. from the points to the tree stem centre
                vector_point_tree_x_list    = x_points_list - Tree_centre_x;
                vector_point_tree_y_list    = y_points_list - Tree_centre_y;

                [x_orth_i_list, y_orth_i_list, delta_orth_list, ~]  = Vector_Circle_Intersection(x_points_list, y_points_list, vector_point_tree_x_list, vector_point_tree_y_list, Tree_centre_x, Tree_centre_y, Tree_radius);

                % These intersection points replace the incorrect values
                x_i_list(missed_beams) = x_orth_i_list(missed_beams);
                y_i_list(missed_beams) = y_orth_i_list(missed_beams);

                % If the points were duplicated, the intersection points must be duplicated as well
                if Point_Duplication == true
                    % The first half belongs to the original points
                    x_i_first_half_list = x_i_list(1 : number_beams/2);
                    y_i_first_half_list = y_i_list(1 : number_beams/2);

                    % The intersection points are duplicated
                    [x_i_duplicated_cell, y_i_duplicated_cell] = Point_Cloud_Duplication(Tree_centre_x, Tree_centre_y, {x_i_first_half_list}, {y_i_first_half_list}, 1);            
                    x_i_list = x_i_duplicated_cell{1};
                    y_i_list = y_i_duplicated_cell{1};

                    % The distances are the same
                    delta_orth_list(number_beams/2 + 1 : end) = delta_orth_list(1 : number_beams/2);
                end

                x_int_points_cell{s} = x_i_list;
                y_int_points_cell{s} = y_i_list;

                %--% The points are mirrored  %--%
                % It is ensured that all points lie on the same side of the circle to avoid the asymmetrical distances affecting the fit by mirroring them
                % The points are moved from inside to outside
                % When moved inside the average distance is lower, and thus the area of non-zero likelihood is greater as well
                outlier_points = delta_orth_list < 0;

                x_points_list(outlier_points) = 2*x_orth_i_list(outlier_points) - x_points_list(outlier_points);
                y_points_list(outlier_points) = 2*y_orth_i_list(outlier_points) - y_points_list(outlier_points);    

                x_points_cell{s} = x_points_list;
                y_points_cell{s} = y_points_list;
            end
        end
        
    %% The expected distance of each point %%
        Expected_distance_cell = cell(1, number_scanners);
        
        for s = 1 : number_scanners
            % This scanner's location and hits
            scanner_x       = Scanner_loc_cell{s}(1);
            scanner_y       = Scanner_loc_cell{s}(2);

            x_points_list   = x_points_cell{s};
            y_points_list   = y_points_cell{s};

            % The uncertainty of these hits
            sigma_radial_list   = sigma_radial_cell{s};
            sigma_prop_list 	= sigma_range_cell{s};

            % The width of these beams at their hit locations
            laser_radius_list   = laser_radius_cell{s};            
            
            %--% The expected distance of each point %--%
            number_beams            = length(x_points_list);
            
            Expected_distance_list  = NaN(1, number_beams);
            
            for b = 1 : number_beams
                % The coordinates and properties of this point
                x_point = x_points_list(b);
                y_point = y_points_list(b);

                sigma_rad       = sigma_radial_list(b);
                sigma_prop      = sigma_prop_list(b);
                laser_radius    = laser_radius_list(b);
                
                % The vector from the scanner to the point
                vector_x    = x_point - scanner_x;
                vector_y    = y_point - scanner_y;
                
                alpha       = pi/2 - atan2(vector_y, vector_x);
                
                % The circle centroid is moved w.r.t. the point
                circle_x_t  = Tree_centre_x - x_point;
                circle_y_t  = Tree_centre_y - y_point;
                
                % And rotated to the (r, p) frame
                circle_r    = cos(alpha) * circle_x_t - sin(alpha) * circle_y_t;
                circle_p    = sin(alpha) * circle_x_t + cos(alpha) * circle_y_t;
                
                % The list of coordinates into which the problem is discretised
                r_list = linspace(-laser_radius, laser_radius, number_steps);
                p_list = linspace(-number_STD * sigma_prop, number_STD * sigma_prop, number_steps);
                
                dr = 2*laser_radius / (number_steps - 1);
                dp = 2*number_STD * sigma_prop / (number_steps - 1);
                
                % If the uncertainty in one direction is zero, the problem becomes one-dimensional
                if sigma_rad == 0
                    % One-dimensional in p
                    L_list = 1/(sqrt(2*pi) * sigma_prop) * exp(-1/2 * (p_list/sigma_prop).^2);
                    
                    D_list = abs(p_list - circle_p - Tree_radius);
                    
                    % The expected_distance
                    Expected_distance = sum(L_list .* D_list.^distance_moment) * dp;                % Note that the distance is to the power of the desired statistical moment
                    
                elseif sigma_prop == 0
                    % One-dimensional in r
                    Z       = erf(laser_radius / (sigma_rad * sqrt(2)));
                    L_list  = 1/(sqrt(2*pi) * sigma_rad * Z) * exp(-1/2 * (r_list/sigma_rad).^2);
                    
                    D_list  = abs(r_list - circle_r - Tree_radius);
                    
                    % The expected distance
                    Expected_distance = sum(L_list .* D_list.^distance_moment) * dr;                % Note that the distance is to the power of the desired statistical moment
                else
                    [r_matrix, p_matrix] = meshgrid(r_list, p_list);

                    % The likelihood of these coordinates
                    Z           = erf(laser_radius / (sigma_rad * sqrt(2)));
                    L_matrix    = 1/(2*pi * sigma_rad * sigma_prop * Z) * exp(-1/2 * ((r_matrix / sigma_rad).^2 + (p_matrix / sigma_prop).^2));

                    % The distances to the circle of these coordinates
                    D_matrix    = sqrt((r_matrix - circle_r).^2 + (p_matrix - circle_p).^2) - Tree_radius;
                    D_matrix    = abs(D_matrix);

                    % The expected distance
                    Expected_distance_matrix    = L_matrix .* D_matrix.^ distance_moment;           % Note that the distance is to the power of the desired statistical moment
                    Expected_distance           = sum(Expected_distance_matrix, 'all') * dr * dp;
                end
                
                Expected_distance_list(b)   = Expected_distance;
            end
            
            Expected_distance_cell{s} = Expected_distance_list;
        end
        
    %% The average expected distance %%
        Expected_distance_list_total    = horzcat(Expected_distance_cell{:});
        Average_expected_distance       = mean(Expected_distance_list_total);

        if Print == true
            fprintf('The average expected distance (%g moment) is %g mm^%g \n', distance_moment, Average_expected_distance * 1e3 ^ distance_moment, distance_moment);
        end
        
   %% Plot %%
        if Plot == true
            %--% The total data %--%
            x_points_list_total     = vertcat(x_points_cell{:});
            y_points_list_total     = vertcat(y_points_cell{:});   
            
            Expected_distance_list_total = Expected_distance_list_total * 1e3 ^ distance_moment;  % m^dm to mm^dm
                        
            %--% Colourmap for the likelihood values %--%
            n_colours   = 1e3;
            purple_cmap = cbrewer('seq', 'Purples', 2*n_colours);
            purple_cmap(1 : n_colours, :) = [];       % The first half of the colours is removed s.t. white dots aren't present

            %--% Coordinates of the tree stem %--%
            theta_list          = linspace(0, 2*pi, number_steps);
            circle_x_list       = Tree_radius * sin(theta_list);
            circle_y_list       = Tree_radius * cos(theta_list);
            
            x_tree_stem_list    = circle_x_list + Tree_centre_x;
            y_tree_stem_list    = circle_y_list + Tree_centre_y;
            
            %--% Plot %--%
            figure(1)

            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
            set(gcf, 'color', [1, 1, 1])
        
            hold on
            grid on

            % Tree stem
            plot(x_tree_stem_list, y_tree_stem_list, 'color', 'k', 'LineWidth', 2, 'DisplayName', 'Tree stem');

            % Point cloud
            if Point_Mirroring == true
                Point_Cloud_String = sprintf('Mirrored point cloud');
            else
                Point_Cloud_String = sprintf('Point cloud');              
            end

            scatter(x_points_list_total, y_points_list_total, [], Expected_distance_list_total, 'filled', 'DisplayName', Point_Cloud_String); 

            % Colourmap
            colormap(purple_cmap);

            cb = colorbar;
            caxis manual
            shading interp
            caxis([0, max(Expected_distance_list_total)])
            colortitlehandle = get(cb, 'Title');
            titlestring = sprintf('Expected distance [mm^%g]', distance_moment);
            set(colortitlehandle, 'String', titlestring);
            cb.FontSize = 15;

            % Legend
            legend('show', 'Location', 'Northoutside');

            % Axis labels
            xlabel('x [m]');
            ylabel('y [m]');

            % Axis limits
            x_lim = Tree_centre_x + 1.5 * Tree_radius * [-1, 1];
            y_lim = Tree_centre_y + 1.5 * Tree_radius * [-1, 1];

            xlim(x_lim);
            ylim(y_lim);

            % The aspect ratio
            AR = (max(y_lim) - min(y_lim)) / (max(x_lim) - min(x_lim));
            pbaspect([1, AR, 1])

            % Axis looks
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
            
            hold off
            
            disp('The script will finish and the plots will close when a key is pressed.')
            pause()

            close(1);
        end

end

