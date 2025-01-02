% This script computes the total likelihood of a circle, given a point cloud
% Note that it expects any bias in the data to already be compensated for

function [Avg_L_N_total, L_N_total_cell, Avg_L_N_max, L_N_max_cell, Avg_L_N_min, L_N_min_cell] = Total_Circle_Likelihood(x_points_cell, y_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, circle_x, circle_y, circle_radius, Scanner_loc_cell, Diagnostics, Point_Mirroring)
   
    %% Inputs %%
        number_theta_steps = 1e4;      % The level of discretisation of the trapezoidal integral

    %% The total likelihood of each point %%
        % The likelihood is determined over the full circle
        theta_list = linspace(0, 2*pi, number_theta_steps);
    
        number_scanners = length(Scanner_loc_cell);
        
        L_N_total_cell  = cell(1, number_scanners);
        L_N_max_cell    = cell(1, number_scanners);
        L_N_min_cell    = cell(1, number_scanners);
                
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
            
            %--% The normalised likelihood of each point %--%
            number_beams    = length(x_points_list);
            
            L_N_total_list  = zeros(1, number_beams);
            L_N_max_list    = zeros(1, number_beams);
            L_N_min_list    = zeros(1, number_beams);

            for b = 1 : number_beams
                % This point's location
                x_point = x_points_list(b);
                y_point = y_points_list(b);
                
                % This point's uncertainty
                sigma_rad       = sigma_radial_list(b);
                sigma_prop      = sigma_prop_list(b);
                laser_radius    = laser_radius_list(b);
                
                L_N_steps = Circle_Likelihood(theta_list, x_point, y_point, scanner_x, scanner_y, circle_x, circle_y, circle_radius, sigma_rad, sigma_prop, laser_radius);
                                
                % The trapezoidal integral value
                L_N_total           = 2*pi/number_theta_steps * sum(L_N_steps);
                L_N_total_list(b)   = L_N_total;
                
                % The maximum
                L_N_max             = max(L_N_steps);
                L_N_max_list(b)     = L_N_max;
                
                % The minimum
                L_N_min             = min(L_N_steps);
                L_N_min_list(b)     = L_N_min;
            end
            
            % The values for this scanner are appended
            L_N_total_cell{s}   = L_N_total_list;
            L_N_max_cell{s}     = L_N_max_list;
            L_N_min_cell{s}     = L_N_min_list;
        end
        
    %% The average likelihood values %%
        Likelihood_N_list_total     = horzcat(L_N_total_cell{:});
        Avg_L_N_total               = mean(Likelihood_N_list_total);
        
        Likelihood_N_list_max       = horzcat(L_N_max_cell{:});
        Avg_L_N_max                 = mean(Likelihood_N_list_max);
        
        Likelihood_N_list_min       = horzcat(L_N_min_cell{:});
        Avg_L_N_min                 = mean(Likelihood_N_list_min);
        
    %% Figure showing the likelihood of each point %%
    if Diagnostics == true
        % The total data
        x_points_list_total     = vertcat(x_points_cell{:});
        y_points_list_total     = vertcat(y_points_cell{:});
        
        Likelihood_Data_Cell    = {Likelihood_N_list_total, Likelihood_N_list_max, Likelihood_N_list_min};
        Likelihood_Names        = {'Total', 'Max', 'Min'};
        
        % Colourmap for the likelihood values
        n_colours   = 1e2;
        purple_cmap = cbrewer('seq', 'Purples', 4*n_colours);
        purple_cmap(1 : n_colours, :) = [];      % The first quarter of the colours is removed s.t. white dots aren't present

        % Coordinates of the tree stem
        x_tree_stem_list = circle_radius .* sin(theta_list) + circle_x;
        y_tree_stem_list = circle_radius .* cos(theta_list) + circle_y;

        %--% Plot %--%
        figure(1)

        % Set the size and white background color
        set(gcf, 'Units', 'Normalized', 'Position', [0 0 0.8 0.8])
        set(gcf, 'color', [1, 1, 1])
        
        for p = 1 : 3
            subplot(1, 3, p);
            
            % This subplot's data
            Likelihood_Data = Likelihood_Data_Cell{p};
            Likelihood_Type = Likelihood_Names{p};
            
            title(Likelihood_Type);
            
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

            scatter(x_points_list_total, y_points_list_total, [], Likelihood_Data, 'filled', 'DisplayName', Point_Cloud_String);

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

            hold off
        end
            
        disp('The script will continue (and close the plot) if you press a key on the keyboard')
        pause
        
        close(1)            
    end
end