% This script determines the maximum normalised likelihood of a circle fit, given a point
% Additionally, the sol_oob flag is true if the solution was not found within the solution interval

function [Avg_L_N_max, L_N_max_cell] = Max_Likelihood_Finder(x_points_cell, y_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, circle_x, circle_y, circle_radius, Scanner_loc_cell, x_int_points_cell, y_int_points_cell, Print, Diagnostics, Point_Mirroring)
   
    %% Inputs %%
        Interval_factor     = 1;        % The factor by which the interval width is multiplied, as a precaution    
        L_N_precision       = 0.01;     % The required precision of the normalised likelihood
        number_theta_steps  = 1e6;      % The level of discretisation to determine the interval width

    %% The maximum likelihood is computed for each point separately %%
        number_scanners     = length(Scanner_loc_cell);
        
        L_N_max_cell        = cell(1, number_scanners);
        theta_optimal_cell  = cell(1, number_scanners);
                
        for s = 1 : number_scanners
            % This scanner's location and hits
            scanner_x       = Scanner_loc_cell{s}(1);
            scanner_y       = Scanner_loc_cell{s}(2);

            x_points_list   = x_points_cell{s};
            y_points_list   = y_points_cell{s};

            % The uncertainty of these hits
            sigma_radial_list   = sigma_radial_cell{s};
            sigma_range_list 	= sigma_range_cell{s};

            % The width of these beams at their hit locations
            laser_radius_list   = laser_radius_cell{s};
            
            % The intersection points
            x_int_points_list   = x_int_points_cell{s};
            y_int_points_list   = y_int_points_cell{s};

            %--% The interval within which the maximum likelihood is found %--%
            number_beams    = length(x_points_list);

            dtheta          = 2*pi/number_theta_steps;
            theta_list      = 0 : dtheta : 2*pi;

            b = 0;
            interval_width = [];     

            while isempty(interval_width) && b < number_beams
                b = b + 1;

                x_point     = x_points_list(b);
                y_point     = y_points_list(b);

                sigma_radial    = sigma_radial_list(b);
                sigma_range     = sigma_range_list(b);
                laser_radius    = laser_radius_list(b);

                % The likelihood values
                L_N_list = Circle_Likelihood(theta_list, x_point, y_point, scanner_x, scanner_y, circle_x, circle_y, circle_radius, sigma_radial, sigma_range, laser_radius);

                % The interval is the width of likelihoods greater than the precision requirement, multiplied by a margin
                first_ind   = find(L_N_list > L_N_precision, 1);
                theta_first = theta_list(first_ind);
                last_ind    = find(L_N_list > L_N_precision, 1, 'last');
                theta_last  = theta_list(last_ind);

                interval_width = (theta_last - theta_first) * Interval_factor;
            end

            % If a zero-likelihood is found for all beams, the closest intersection is used instead
            if isempty(interval_width)
                interval_width = 0;
            end            
        
            %--% The initial guesses for the Golden-section search %--%
            % The angle for the initial guess is the angle from the circle centre to the intersection point
            delta_x_list = x_int_points_list - circle_x;
            delta_y_list = y_int_points_list - circle_y;

            theta_intersect_list = pi/2 - atan2(delta_y_list, delta_x_list);

            %--% Golden-section search to find the highest likelihood %--%
            % The golden ratio
            psi = (1 + sqrt(5))/2;

            L_N_max_list        = zeros(1, number_beams);
            theta_optimal_list  = zeros(1, number_beams);
            number_oob          = 0;     % The number of times the optimum is outside the search interval is saved

            for b = 1 : number_beams
                % This point
                x_point = x_points_list(b);
                y_point = y_points_list(b);
                
                % The closest intersection points are used to establish the search interval
                theta_2_init    = theta_intersect_list(b);     % One of the central points is given the intersect value

                theta_1_init    = theta_2_init - (2 - psi) * interval_width;
                theta_4_init    = theta_1_init + interval_width;

                theta_1 = theta_1_init;
                theta_4 = theta_4_init;

                % Convergence continues until the change in likelihood meets the desired precision threshold
                convergence = false;

                while convergence == false
                    % The points at the golden ratios (from x1 and x4)
                    theta_2 = theta_1 + (2 - psi) * (theta_4 - theta_1);
                    theta_3 = theta_1 + (psi - 1) * (theta_4 - theta_1);

                    % Likelihood values at theta 1 to 4
                    theta_GS_list = [theta_1, theta_2, theta_3, theta_4];

                    L_N_GS_list = Circle_Likelihood(theta_GS_list, x_point, y_point, scanner_x, scanner_y, circle_x, circle_y, circle_radius, sigma_radial, sigma_range, laser_radius);

                    % The point of maximum likelihood
                    L_N_max         = max(L_N_GS_list);
                    L_N_max_point   = find(L_N_GS_list == L_N_max, 1);

                    % Convergence check    
                    L_N_min = min(L_N_GS_list);     

                    if L_N_max - L_N_min < L_N_precision    % If the amplitude is smaller than the threshold, precision must exceed it
                        convergence = true;

                    elseif theta_4 - theta_1 < 1e-6         % Alternatively if the gap between angles becomes so small, the script continues
                        convergence = true;                 % This can happen if the highest normalised likelihood is lower than the required precision

                    % If convergence is not met, the points are updated
                    else                    
                        if L_N_max_point <= 2
                            theta_4 = theta_3;
                        else
                            theta_1 = theta_2;
                        end
                    end
                end

                % The maximum likelihood is appended
                L_N_max_list(b) = L_N_max;

                % Check if the bounds were reached
                theta_optimal           = theta_GS_list(L_N_max_point);
                theta_optimal_list(b)   = theta_optimal;

                % If the likelihood is near-zero, the outer bounds are more likely to be hit due to the poor solution space and hitting them is less relevant
                % Additionally, this is only relevant if the interval width was nonzero        
                if L_N_max > L_N_precision && interval_width > 0
                    if theta_optimal == theta_1_init || theta_optimal == theta_4_init
                        number_oob = number_oob + 1;
                    end
                end      
                
                %--% Show diagnostic information to check the quality of the results %--%
                if Diagnostics == true
                    % Normalised likelihood in both directions
                    [r, p]  = Coordinate_Translator(theta_optimal, circle_radius, circle_x, circle_y, x_point, y_point, [scanner_x, scanner_y]);

                    L_N_radial      = Radial_Likelihood(r, sigma_radial, laser_radius);
                    L_N_propagation = Propagation_Likelihood(p, sigma_range);
                        
                    disp('------------------------------')
                    fprintf('The normalised maximum likelihood of this point is: %g \n', L_N_max);
                    fprintf('   Norm. likelihood radial: %g \n', L_N_radial);
                    fprintf('   Norm. likelihood propagation: %g \n', L_N_propagation);
                    disp('------------------------------')                
                    fprintf('The initial estimate of the angle (closest intersect) is: %g degrees \n', rad2deg(theta_2_init));
                    fprintf('The final estimate of the angle is: %g degrees \n', rad2deg(theta_optimal));
                    disp('------------------------------')                
                    fprintf('The standard deviation in the radial direction is: %g mm \n', sigma_radial * 1e3);
                    fprintf('The standard deviation in the propagation direction is: %g mm \n', sigma_range * 1e3);                
                    disp('------------------------------')                

                    % Locations of the initial and final estimates
                    x_init = circle_radius * sin(theta_2_init) + circle_x;
                    y_init = circle_radius * cos(theta_2_init) + circle_y;

                    x_optimal = circle_radius * sin(theta_optimal) + circle_x;
                    y_optimal = circle_radius * cos(theta_optimal) + circle_y;

                    % Plot
                    figure(1)

                    % Set the size and white background color
                    set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
                    set(gcf, 'color', [1, 1, 1])

                    hold on
                    grid on

                    % Tree stem
                    Tree_stem = viscircles([circle_x, circle_y], circle_radius, 'color', 'k', 'LineWidth', 2);

                    x_t     = Tree_stem.Children(1).XData;
                    x_t     = x_t(~isnan(x_t));
                    y_t     = Tree_stem.Children(1).YData;
                    y_t     = y_t(~isnan(y_t));

                    Tree    = fill(x_t, y_t, 'k', 'FaceAlpha', 0.10);
                    Tree.DisplayName = 'Tree stem';

                    % Laser beam
                    plot([scanner_x, x_point], [scanner_y, y_point], 'color', 'r', 'LineWidth', 2, 'DisplayName', 'Laser beam');

                    % Initial estimate
                    scatter(x_init, y_init, 'DisplayName', 'Initial estimate', 'MarkerFaceColor', 'b');

                    % Final estimate
                    scatter(x_optimal, y_optimal, 'DisplayName', 'Final estimate', 'MarkerFaceColor', 'g');

                    % Legend
                    legend('show', 'Location', 'Northoutside');

                    % Axis labels
                    xlabel('x [m]');
                    ylabel('y [m]');

                    % Axis limits
                    x_lim = circle_x + 1.5 * circle_radius * [-1, 1];
                    y_lim = circle_y + 1.5 * circle_radius * [-1, 1];

                    xlim(x_lim);
                    ylim(y_lim);

                    % The aspect ratio
                    AR = (max(y_lim) - min(y_lim)) / (max(x_lim) - min(x_lim));
                    pbaspect([1, AR, 1])

                    % Axis looks
                    set(gca, 'FontSize', 15);
                    set(gca, 'LineWidth', 2);

                    disp('The script will continue (and close the plot) if you press a key on the keyboard')
                    pause

                    hold off
                    close(1)
                end                
            end
            
            % The results for this scanner
            L_N_max_cell{s}         = L_N_max_list;
            theta_optimal_cell{s}   = theta_optimal_list;
            
            % The percentage of beams for which the solution lay outside the interval
            percentage_oob = number_oob / number_beams * 100;

            if Print == true && percentage_oob > 0
                fprintf('   The optimum was not found within the interval for %.2g%% of the beams \n', percentage_oob);
                disp('  If the percentage is high or this message is printed often, consider increasing the factor by which the GSS interval is multiplied')
            end
        end
        
    %% The average maximum likelihood %%
        Likelihood_N_list_total = horzcat(L_N_max_cell{:});
        Avg_L_N_max             = mean(Likelihood_N_list_total);
        
    %% Figure showing the likelihood of each point %%
    if Diagnostics == true
        % The total data
        theta_opt_list_total    = horzcat(theta_optimal_cell{:});
        x_points_list_total     = vertcat(x_points_cell{:});
        y_points_list_total     = vertcat(y_points_cell{:});
        
        % Colourmap for the likelihood values
        n_colours   = 1e3;
        purple_cmap = cbrewer('seq', 'Purples', 2*n_colours);
        purple_cmap(1 : n_colours, :) = [];      % The first half of the colours is removed s.t. white dots aren't present

        % Coordinates of the optimal points
        x_opt_list = circle_radius .* sin(theta_opt_list_total) + circle_x;
        y_opt_list = circle_radius .* cos(theta_opt_list_total) + circle_y;

        % Coordinates of the tree stem
        x_tree_stem_list = circle_radius .* sin(theta_list) + circle_x;
        y_tree_stem_list = circle_radius .* cos(theta_list) + circle_y;

        % Plot
        figure(2)

        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
        set(gcf, 'color', [1, 1, 1])

        for p = 1 : 2
            subplot(1, 2, p)

            hold on
            grid on

            % Tree stem
            plot(x_tree_stem_list, y_tree_stem_list, 'color', 'k', 'LineWidth', 2, 'DisplayName', 'Tree stem');

            % Point cloud
            if p == 1
                if Point_Mirroring == true
                    Point_Cloud_String = sprintf('Mirrored point cloud');
                else
                    Point_Cloud_String = sprintf('Point cloud');              
                end
                
                scatter(x_points_list_total, y_points_list_total, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', 'DisplayName', Point_Cloud_String, 'MarkerFaceAlpha', 1);
            end

            % The optimal points
            if p == 2
                scatter(x_opt_list, y_opt_list, [], Likelihood_N_list_total, 'filled', 'DisplayName', 'Optimal points');

                % Colourmap
                colormap(purple_cmap);

                cb = colorbar;
                caxis manual
                shading interp
                caxis([0, 1])
                colortitlehandle = get(cb, 'Title');
                titlestring = 'L_{N} [-]';
                set(colortitlehandle, 'String', titlestring);
                cb.FontSize = 15;
            end

            % Legend
            legend('show', 'Location', 'Northoutside');

            % Axis labels
            xlabel('x [m]');
            ylabel('y [m]');

            % Axis limits
            x_lim = circle_x + 1.5 * circle_radius * [-1, 1];
            y_lim = circle_y + 1.5 * circle_radius * [-1, 1];

            xlim(x_lim);
            ylim(y_lim);

            % The aspect ratio
            AR = (max(y_lim) - min(y_lim)) / (max(x_lim) - min(x_lim));
            pbaspect([1, AR, 1])

            % Axis looks
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
        end

        disp('The script will continue (and close the plot) if you press a key on the keyboard')
        pause

        hold off
        close(2)            
    end 
end