% This script filters the point cloud for outliers by using a simple standard deviation threshold
% Each scanner is assessed individually, as the uncertainty of different scanners can vary

% Iteration can be set to true or false to indicate whether or not the geometry and outliers are allowed to change iteratively until convergence
% If it is set to false, the circle geometry must be given in the form [x, y, r]

function [x_f_points_cell, y_f_points_cell, number_beams_list] = Point_Cloud_Filtering(x_points_cell, y_points_cell, number_STD, Iteration, Input_Geometry, Weighting, Scanner_loc_cell, beam_divergence, sigma_range_device, range_bias, max_incidence_angle, Print, Plot)

    %% Inputs %%
        convergence_threshold   = 001;  % [%] If the change in geometry falls below this value, iteration ends
        max_iter                = 1e2;  % Alternatively, the maximum number of iterations also limits it
        
    %% Check whether or not the geometry was given %%
        if isempty(Input_Geometry)
            Geometry_Given  = false;
        else
            Geometry_Given  = true;
            circle_x        = Input_Geometry(1);
            circle_y        = Input_Geometry(2);
            circle_r        = Input_Geometry(3);
        end
        
        if Iteration == false && Geometry_Given == false
            disp('Error: Filtering can not be performed without a given geometry or iteration')
        end
        
    %% Point cloud filtering %%
        number_scanners     = length(x_points_cell);

        x_f_points_cell     = cell(1, number_scanners);
        y_f_points_cell     = cell(1, number_scanners);
        
        number_beams_list   = zeros(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's point cloud and location
            scanner_loc     = Scanner_loc_cell{s};
            x_points_list   = x_points_cell{s};
            y_points_list   = y_points_cell{s};
            
            %--% The given geometry is used to determine the outliers %--%
            if      Iteration == false
                [x_f_points_list, y_f_points_list, radial_outliers, angular_outliers, outliers, delta_theta_list, theta_filtering_thresholds, delta_radius_list, radius_filtering_thresholds] = Filtering_Function(x_points_list, y_points_list, circle_x, circle_y, circle_r, number_STD);
                
            %--% Iteration until the geometry (and outliers) have converged %--%
            elseif  Iteration == true
            
                % If the geometry is not given, an initial guess is made
                if      Geometry_Given == false
                    circle_x_old    = mean(x_points_list);
                    circle_y_old    = mean(y_points_list);
                    circle_r_old    = mean(sqrt((x_points_list - circle_x_old).^2 + (y_points_list - circle_y_old).^2));
                    
                % Otherwise it is simply used as the initial guess
                elseif  Geometry_Given == true
                    circle_x_old    = circle_x;
                    circle_y_old    = circle_y;
                    circle_r_old    = circle_r;
                end

                convergence     = false;
                iter            = 0;

                while convergence == false && iter < max_iter
                    iter            = iter + 1;

                    % Point cloud filtering
                    [x_f_points_list, y_f_points_list, radial_outliers, angular_outliers, outliers, delta_theta_list, theta_filtering_thresholds, delta_radius_list, radius_filtering_thresholds] = Filtering_Function(x_points_list, y_points_list, circle_x_old, circle_y_old, circle_r_old, number_STD);                    
                    
                    % New geometry is fitted based on this filtered point cloud %
                    Method              = 'LS';
                    Point_Duplication   = false;
                    [circle_x, circle_y, circle_r] = Least_Squares_Circle_Fitting({x_f_points_list}, {y_f_points_list}, Method, Weighting, {scanner_loc}, beam_divergence, sigma_range_device, range_bias, max_incidence_angle, Point_Duplication, []);

                    % For convergence only the radius is checked, as the x and y coordinates can be near-zero
                    change = (circle_r - circle_r_old) / circle_r_old * 100;

                    if abs(change) < convergence_threshold
                        convergence = true;

                    else
                        circle_x_old = circle_x;
                        circle_y_old = circle_y;
                        circle_r_old = circle_r;
                    end
                end
            end
            
            % The filtered points are appended for this scanner
            x_f_points_cell{s} = x_f_points_list;
            y_f_points_cell{s} = y_f_points_list;
            
            % The new number of beams
            number_beams            = length(x_f_points_list);
            number_beams_list(s)    = number_beams;
            
            % The percentage of outliers is printed
            if Print == true
                percentage_outliers = sum(outliers) / length(x_points_list) * 100;
                fprintf('The percentage of outliers of scanner %g was %.2g%% \n', s, percentage_outliers);
            end
            
            % Plot showing the filtering results
            if Plot == true
                number_beams_original = length(x_points_list);
                
                figure(313)
                
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
                set(gcf, 'color', [1, 1, 1])
                
                %--% Point cloud %--%
                subplot(1, 3, 1)
                
                hold on
                grid on
                
                scatter(x_points_list, y_points_list, 'MarkerFaceColor', 'b', 'DisplayName', 'Original point cloud');
                scatter(x_f_points_list, y_f_points_list, 'MarkerFaceColor', 'r', 'DisplayName', 'Filtered point cloud');
                
                % Legend
                legend('show', 'Location', 'Northoutside');
                
                % Axis labels
                xlabel('x [m]');
                ylabel('y [m]');
                
                % Axis limits
                x_lim = [min(x_points_list) - 0.1, max(x_points_list) + 0.1];
                y_lim = [min(y_points_list) - 0.1, max(y_points_list) + 0.1];
                
                xlim(x_lim);
                ylim(y_lim);

                % The aspect ratio
                AR = (max(y_lim) - min(y_lim)) / (max(x_lim) - min(x_lim));
                pbaspect([1, AR, 1])
                            
                % Axis looks
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
                hold off
                
                %--% Radius %--%
                subplot(1, 3, 2)
                hold on
                grid on
                
                scatter(1 : number_beams_original, delta_radius_list, 'DisplayName', 'Filtered point cloud', 'MarkerFaceColor', 'r');                
                scatter(find(radial_outliers == true), delta_radius_list(radial_outliers), 'DisplayName', 'Original point cloud', 'MarkerFaceColor', 'b');
                
                plot([0, number_beams_original], [min(radius_filtering_thresholds), min(radius_filtering_thresholds)], 'DisplayName', 'Lower threshold', 'LineWidth', 2, 'color', 'm');
                plot([0, number_beams_original], [max(radius_filtering_thresholds), max(radius_filtering_thresholds)], 'DisplayName', 'Upper threshold', 'LineWidth', 2, 'color', 'c');
                
                % Legend
                legend('show', 'Location', 'Northoutside');
                
                % Axis labels
                xlabel('Beam [-]');
                ylabel('\Delta R [m]');
                
                % Axis looks
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
                hold off
                
                xlim([0, length(delta_radius_list)]);
                
                %--% Theta %--%
                subplot(1, 3, 3)
                hold on
                grid on
                
                scatter(1 : number_beams_original, rad2deg(delta_theta_list), 'DisplayName', 'Filtered point cloud', 'MarkerFaceColor', 'r');
                scatter(find(angular_outliers == true), rad2deg(delta_theta_list(angular_outliers)), 'DisplayName', 'Original point cloud', 'MarkerFaceColor', 'b');

                plot([0, number_beams_original], rad2deg([max(theta_filtering_thresholds), max(theta_filtering_thresholds)]), 'DisplayName', 'Upper threshold', 'LineWidth', 2, 'color', 'c');
                
                % Legend
                legend('show', 'Location', 'Northoutside');
                
                % Axis labels
                xlabel('Beam [-]');
                ylabel('\Delta \Theta [deg]');
                
                % Axis looks
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
                hold off
                
                xlim([0, length(delta_theta_list)]);
                
                disp('The script will continue (and close the plot) if you press a key on the keyboard')
                pause

                hold off
                close(313)
            end
        end
        
    %% Filtering function %%
    function [x_f_points_list, y_f_points_list, radial_outliers, angular_outliers, outliers, delta_theta_list, theta_filtering_thresholds, delta_radius_list, radius_filtering_thresholds] = Filtering_Function(x_points_list, y_points_list, circle_x, circle_y, circle_r, number_STD)
        % The centered point cloud
        delta_x_list    = x_points_list - circle_x;
        delta_y_list    = y_points_list - circle_y;

        %--% Angular outliers %--%
        % The angles as viewed from the centre of the circle                        
        theta_list          = atan(delta_x_list ./ delta_y_list);
        delta_theta_list    = theta_list - mean(theta_list);

        % The theta filtering thresholds
        theta_STD = std(delta_theta_list);

        theta_filtering_thresholds = mean(delta_theta_list) + number_STD * theta_STD * [-1, 1];

        % The angular outliers
        angular_outliers_low    = delta_theta_list < min(theta_filtering_thresholds);
        angular_outliers_high   = delta_theta_list > max(theta_filtering_thresholds);

        angular_outliers        = boolean(angular_outliers_low + angular_outliers_high);

        %--% Radial outliers %--%
        % The resulting radii
        radius_list         = sqrt(delta_x_list.^2 + delta_y_list.^2);
        delta_radius_list   = radius_list - circle_r;   % The radius of the circle is subtracted 

        % The radius filtering thresholds
        radius_STD = std(delta_radius_list);

        radius_filtering_thresholds = mean(delta_radius_list) + number_STD * radius_STD * [-1, 1];

        % The radial outliers
        radius_outliers_low     = delta_radius_list < min(radius_filtering_thresholds);
        radius_outliers_high    = delta_radius_list > max(radius_filtering_thresholds);

        radial_outliers         = boolean(radius_outliers_low + radius_outliers_high);

        %--% Filtering of all outliers %--%
        outliers                = boolean(angular_outliers + radial_outliers);

        % The filtered point cloud is returned
        x_f_points_list = x_points_list(~outliers);
        y_f_points_list = y_points_list(~outliers);
    end
end