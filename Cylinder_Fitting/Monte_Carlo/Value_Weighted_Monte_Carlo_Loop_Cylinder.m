% The value of each point is used to weigh their contribution to cylinder fitting, which is the ratio of its information to said information's uncertainty
% These metrics are first determined using leave-one-out unweighted Monte Carlo loops

function [Cylinder_Geometry, Gaussian_Mixture_Models, Monte_Carlo_Loop, Fitting_Parameters] = Value_Weighted_Monte_Carlo_Loop_Cylinder(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Monte_Carlo_Inputs, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions)

    %% Inputs %%
        
        % Point weighting
        Point_Weighting         = Fitting_Parameters.Point_Weighting;

        % Point cloud
        point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;
        number_points_list      = Point_Cloud_Coord.number_points_list;

        % Scanners
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;

        % Scanner parameters
        range_bias              = Scanner_Parameters.range_bias;

        % Statistical values
        Confidence_interval     = Statistical_Values.Confidence_interval;

        % Monte Carlo
        max_MC_length           = Monte_Carlo_Inputs.max_MC_length;

        % Outputs
        Plot                    = Output_Decisions.Plot;
        
    %% Uniform weights %%
        if strcmp(Point_Weighting, 'None')
            % In case no weighting is desired, it is equal for each point
            number_points       = sum(number_points_list);
            point_weights_list  = ones(number_points, 1);    
        end
        
    %% Information %%        
        if strcmp(Point_Weighting, 'Info') ||strcmp(Point_Weighting, 'Value')
            % The information of each point is seen as the extent to which the geometry deviates when it is left out
            Output_Decisions.Progress_Update    = true;
            [Information_list, ~, ~]            = Information_Cylinder_Point_Cloud(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions);
            point_weights_list                  = Information_list;
        end

    %% Reliability %%
        if strcmp(Point_Weighting, 'Reliability') || strcmp(Point_Weighting, 'Value')
            %--% Reliability %--%            
            % In order to sample the point cloud, its uncertainty must be estimated which is a product of the geometry
            % Outputs are recommended to be off
            LS_Diagnostics                          = false;                % [true, false]
            Sampling_Diagnostics                    = false;                % [true, false]
            
            Output_Decisions_MC                     = Output_Decisions;
            Output_Decisions_MC.Plot                = false;                % [true, false]
            Output_Decisions_MC.Progress_Update     = false;                % [true, false]

            % The initial geometry is determined through least-squares using a matrix of the point cloud    
            LS_Geometry = Least_Squares_Cylinder_Fitting(Point_Cloud_Coord, Scanning_Parameters, Fitting_Parameters, Confidence_interval, LS_Diagnostics, LS_Diagnostics);
    
            cylinder_centre_init    = LS_Geometry.Cylinder_centre;
            cylinder_direction_init = LS_Geometry.Cylinder_direction;
            cylinder_radius_init    = LS_Geometry.Cylinder_radius;
            cylinder_length_init    = LS_Geometry.Cylinder_length;

            % The uncertainty of each point, given the geometry
            [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(cylinder_centre_init, cylinder_direction_init, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
            alpha                       = 1 - Confidence_interval/100;
            Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Sampling_Diagnostics);

            % The uncertainty of each point's information
            Information_cell = cell(1, max_MC_length);

            for i = 1 : max_MC_length
                % Sample the uncertain point cloud
                Point_Cloud_Coord_MC        = Gaussian_Point_Sampling(Point_Cloud_Distributions, Sampling_Diagnostics, cylinder_centre_init, cylinder_radius_init, cylinder_direction_init);

                % The information of this sampled point cloud
                [Information_list, ~, ~]    = Information_Cylinder_Point_Cloud(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord_MC, Output_Decisions_MC);
                Information_cell{i}         = Information_list;

                % Progress
                progress        = floor(i / max_MC_length * 100);
                progress_last   = floor((i - 1) / max_MC_length * 100);

                if progress - progress_last >= 1
                    fprintf('Reliability estimation progress: %i%% \n', progress);
                end
            end

            % The information uncertainty and resulting reliability
            Information_matrix  = horzcat(Information_cell{:});
            sigma_info_list     = std(Information_matrix, [], 2);

            Reliability_list    = mean(sigma_info_list) ./ sigma_info_list;
            point_weights_list  = Reliability_list;

            %--% Plot %--%
            if Plot == true
                % Colour map
                num_colours         = 1e3;
                reliability_cmap    = cbrewer('seq', 'Greens', num_colours);
                reliability_cmap    = max(reliability_cmap, 0);
                reliability_cmap    = min(reliability_cmap, 1);

                reliability_min     = min(Reliability_list);
                reliability_max     = max(Reliability_list);
                cmap_ind            = round((num_colours - 1) * (Reliability_list - reliability_min)/(reliability_max - reliability_min)) + 1;  
    
                % The least-squares cylinder geometry for reference 
                num_points_cyl = 1e3;
                [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius_init, cylinder_length_init, cylinder_centre_init, cylinder_direction_init, num_points_cyl);

                figure(50);
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])     
        
                % The cylinder on the first set of axes
                ax1 = axes;
                grid on
                hold on

                cyl_name = 'Cylinder';
                cyl_surf = surf(ax1, cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', cyl_name);
        
                % Point cloud shaded by the data in the second set of axes        
                ax2 = axes;
                grid on
                hold on
                pt_name             = 'Point cloud';
                point_cloud_matrix  = vertcat(point_cloud_cell{:});
                pt_sc               = scatter3(ax2, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 200, reliability_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', pt_name);
        
                % Colorbar (on second axes)
                colormap(ax2, reliability_cmap);
        
                cb = colorbar(ax2);
                shading interp
                clim([reliability_min, reliability_max])
                ylabel(cb, 'R [-]');
                cb.FontSize = 25;        
        
                % Axes
                xlabel(ax1, 'x [m]')
                ylabel(ax1, 'y [m]')
                zlabel(ax1, 'z [m]')
        
                x_lim = [min(cylinder_coord_x(:)) - cylinder_radius_init, max(cylinder_coord_x(:)) + cylinder_radius_init];
                y_lim = [min(cylinder_coord_y(:)) - cylinder_radius_init, max(cylinder_coord_y(:)) + cylinder_radius_init];
                z_lim = [min(cylinder_coord_z(:)) - cylinder_radius_init, max(cylinder_coord_z(:)) + cylinder_radius_init];
        
                xlim([ax1, ax2], x_lim);
                ylim([ax1, ax2], y_lim);
                zlim([ax1, ax2], z_lim);
        
                ax1.DataAspectRatio = [1, 1, 1];
                ax2.DataAspectRatio = [1, 1, 1];
                view([ax1, ax2], 45, 45)
                
                % Legend
                legend(ax1, [cyl_surf, pt_sc], {cyl_name, pt_name}, 'location', 'northoutside');
        
                set([ax1, ax2], 'FontSize', 15);
                set([ax1, ax2], 'LineWidth', 2);
        
                % Change the size, link the axes and make the second set invisible
                axis_dimensions     = [ax1.Position; ax2.Position];
                axis_starts         = max(axis_dimensions(:, 1:2), [], 1);
                axis_sizes          = min(axis_dimensions(:, 3:4), [], 1);    
                ax1.Position        = [axis_starts, axis_sizes];
                ax2.Position        = [axis_starts, axis_sizes];
                ax2.Visible         = 'off';
                
                linkprop([ax1, ax2], {'View', 'XLim', 'YLim', 'ZLim'});

                hold([ax1, ax2], 'off');        

                % The figure is saved
                figure_name = 'Point_Cloud_Reliability.png';
                export_fig(50, figure_name);

                if strcmp(Point_Weighting, 'Reliability')
                    % Pause message
                    disp('The reliability plot will close and the script will continue upon a key-press.');
                    pause();
        
                    close(50)
                end
            end

        end

    %% Value %%
        if strcmp(Point_Weighting, 'Value')
            %--% Value %--%
            % The value of each point is the product of the information and its reliability
            Value_list          = Information_list .* Reliability_list;
            point_weights_list  = Value_list;

            %--% Plot %--%
            if Plot == true
                % Colour map
                value_cmap  = cbrewer('seq', 'Oranges', num_colours);
                value_cmap  = max(value_cmap, 0);
                value_cmap  = min(value_cmap, 1);

                value_min   = min(Value_list);
                value_max   = max(Value_list);
                cmap_ind    = round((num_colours - 1) * (Value_list - value_min)/(value_max - value_min)) + 1;  

                figure(100);
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])     
        
                % The cylinder on the first set of axes
                ax1 = axes;
                grid on
                hold on

                cyl_name = 'Cylinder';
                cyl_surf = surf(ax1, cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'None', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', cyl_name);
        
                % Point cloud shaded by the data in the second set of axes        
                ax2 = axes;
                grid on
                hold on
                pt_name             = 'Point cloud';
                point_cloud_matrix  = vertcat(point_cloud_cell{:});
                pt_sc               = scatter3(ax2, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 200, value_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', pt_name);
        
                % Colorbar (on second axes)
                colormap(ax2, value_cmap);
        
                cb = colorbar(ax2);
                shading interp
                clim([value_min, value_max])
                ylabel(cb, 'V [-]');
                cb.FontSize = 25;        
        
                % Axes
                xlabel(ax1, 'x [m]')
                ylabel(ax1, 'y [m]')
                zlabel(ax1, 'z [m]')
        
                x_lim = [min(cylinder_coord_x(:)) - unw_cyl_radius, max(cylinder_coord_x(:)) + unw_cyl_radius];
                y_lim = [min(cylinder_coord_y(:)) - unw_cyl_radius, max(cylinder_coord_y(:)) + unw_cyl_radius];
                z_lim = [min(cylinder_coord_z(:)) - unw_cyl_radius, max(cylinder_coord_z(:)) + unw_cyl_radius];
        
                xlim([ax1, ax2], x_lim);
                ylim([ax1, ax2], y_lim);
                zlim([ax1, ax2], z_lim);
        
                ax1.DataAspectRatio = [1, 1, 1];
                ax2.DataAspectRatio = [1, 1, 1];
                view([ax1, ax2], 45, 45)
                
                % Legend
                legend(ax1, [cyl_surf, pt_sc], {cyl_name, pt_name}, 'location', 'northoutside');
        
                set([ax1, ax2], 'FontSize', 15);
                set([ax1, ax2], 'LineWidth', 2);
        
                % Change the size, link the axes and make the second set invisible
                axis_dimensions     = [ax1.Position; ax2.Position];
                axis_starts         = max(axis_dimensions(:, 1:2), [], 1);
                axis_sizes          = min(axis_dimensions(:, 3:4), [], 1);    
                ax1.Position        = [axis_starts, axis_sizes];
                ax2.Position        = [axis_starts, axis_sizes];
                ax2.Visible         = 'off';
                
                linkprop([ax1, ax2], {'View', 'XLim', 'YLim', 'ZLim'});

                hold([ax1, ax2], 'off');        

                % The figure is saved
                figure_name = 'Point_Cloud_Value.png';
                export_fig(100, figure_name);

                % Pause message
                disp('The value plot has been created and saved. The script will continue automatically.');    
                close all
            end
        end

    %% (Un)weighted Monte Carlo loop %%
        % The point weights are appended to the fitting parameters
        Fitting_Parameters.point_weights_list = point_weights_list;

        % (Un)weighted Monte Carlo loop
        [Cylinder_Geometry, Gaussian_Mixture_Models, Monte_Carlo_Loop] = Monte_Carlo_Loop_Cylinder(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Monte_Carlo_Inputs, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions);
            
end