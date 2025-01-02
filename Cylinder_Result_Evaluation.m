% The results of the Monte Carlo loop are shown, in comparison with the true geometry

function Results_Tables = Cylinder_Result_Evaluation(True_Cylinder_Geometry, Point_Cloud_Coord, Optim_Cylinder_Geometry_GMM, Gaussian_Mixture_Models, Monte_Carlo_Loop, Scanning_Parameters, Scanner_Parameters, Fitting_Parameters, Statistical_Values, Output_Decisions)

    %% Inputs %% 
        % The scanners
        number_scanners             = Scanning_Parameters.number_scanners;
        Scanner_loc_cell            = Scanning_Parameters.Scanner_loc_cell;

        % Scanner parameters
        range_bias                  = Scanner_Parameters.range_bias;

        % Point cloud
        point_cloud_cell            = Point_Cloud_Coord.point_cloud_cell;
    
        % The true geometry
        Cylinder_centre             = True_Cylinder_Geometry.Cylinder_centre;
        Cylinder_direction          = True_Cylinder_Geometry.Cylinder_direction;
        Cylinder_radius             = True_Cylinder_Geometry.Cylinder_radius;
        Cylinder_length             = True_Cylinder_Geometry.Cylinder_length;
        
        % Monte Carlo data
        MC_Geometry_Data            = Monte_Carlo_Loop.MC_Geometry_Data;
        number_MC_iterations        = Monte_Carlo_Loop.number_MC_iterations;
        expected_Mahal_dist_matrix  = Monte_Carlo_Loop.expected_Mahal_dist_matrix;
        incidence_angle_matrix      = Monte_Carlo_Loop.incidence_angle_matrix;
        propagation_error_matrix    = Monte_Carlo_Loop.propagation_error_matrix;
        cylinder_distance_matrix    = Monte_Carlo_Loop.cylinder_distance_matrix;

        % Fitting parameters
        distance_moment             = Fitting_Parameters.distance_moment;
        Point_Weighting             = Fitting_Parameters.Point_Weighting;
        point_weights_list          = Fitting_Parameters.point_weights_list;

        % Statistical values
        Confidence_interval         = Statistical_Values.Confidence_interval;

        % Output decisions
        Plot                        = Output_Decisions.Plot;
        Print                       = Output_Decisions.Print;        

    %% Save folder %%
        % A folder is created with the current time as its name to save the results in
        current_time        = datetime('now', 'format', 'yyyyMMdd_HHmmss');                         % Date followed by clock time       
        save_folder_name    = string(current_time);                                                 % Converted to string array that acts as the name of the folder

        mkdir(save_folder_name);
        
        % Inputs and outputs of the cylinder fitting algorithm are saved
        save('Cylinder_Fitting_Inputs', 'Point_Cloud_Coord', 'Scanning_Parameters', 'Scanner_Parameters', 'Fitting_Parameters', 'Statistical_Values', 'True_Cylinder_Geometry');
        save('Cylinder_Fitting_Outputs', 'Optim_Cylinder_Geometry_GMM', 'Gaussian_Mixture_Models', 'Monte_Carlo_Loop');
                
        % Data files and figures created so far are moved to the save folder, if they exist
        data_formats        = {'.mat', '.png', '.fig'};
        number_data_formats = length(data_formats);

        for d = 1 : number_data_formats
            try
                data_format_files = ['*', data_formats{d}];
                movefile(data_format_files, save_folder_name);
            catch
            end
        end
        
    %% Analysed variables %%
        % Cylinder variables
        variable_fields         = {'Cylinder_centre', 'Cylinder_direction', 'Cylinder_radius', 'Cylinder_length', 'Cylinder_height_top', 'Cylinder_height_bot'};
        variable_labels         = {'c_x', 'c_y', 'c_z', 'v_x', 'v_y', 'v_z', 'r', 'l', 'h_t', 'h_l'};           % Note that the centre and direction are multi-dimensional
        number_variables        = length(variable_fields);

        % The relative standard deviation and error are divided by the variables or values given here
        normalisation_fields    = {'Cylinder_radius', 1, 'Cylinder_radius', 'Cylinder_length', 'Cylinder_height_top', 'Cylinder_height_bot'};
        
    %% Optimisation data %%
        % The mean and standard deviations as well as the optimiser steps
        [Fuzzy_MC_Geometry, LS_MC_Geometry, Objective_Steps, Geometry_Steps, Geometry_Gradient_Steps, Step_Sizes] = Cylinder_Optimisation_History(Monte_Carlo_Loop);

        % The infinite cylinder geometry error of each optimiser as well as the number of steps the optimiser took
        [LS_MC_Geometry, Fuzzy_MC_Geometry, Fuzzy_GMM_Geometry, ~, objective_value_list, number_steps_list, number_UU_iter_list] = Monte_Carlo_Loop_Cylinder_Performance(Monte_Carlo_Loop, True_Cylinder_Geometry, LS_MC_Geometry, Fuzzy_MC_Geometry, Optim_Cylinder_Geometry_GMM, variable_fields, normalisation_fields);    

    %% Tabulated estimation errors %%
        % The fitting methods are aggregated into cell arrays for simplicity        
        Fitting_Method_Data_cell    = {LS_MC_Geometry, Fuzzy_MC_Geometry, Fuzzy_GMM_Geometry};
        fitting_method_names        = {'LS (MC)', 'Fuzzy (MC)', 'Fuzzy (GMM)'};
        fitting_method_labels       = {'LS_MC', 'Fuzzy_MC', 'Fuzzy_GMM'};
        number_fit_methods          = length(Fitting_Method_Data_cell);

        for f = 1 : number_fit_methods
            % This fitting method's data
            fitting_method      = fitting_method_names{f};
            fitting_label       = fitting_method_labels{f};
            Fitting_Method_Data = Fitting_Method_Data_cell{f};

            % Data for each variable
            [variable_true_cell, variable_mean_cell, variable_rel_error_cell, variable_STD_cell, variable_rel_STD_cell] = deal(cell(1, number_variables));

            for v = 1 : number_variables
                % This variable's (relativisation) field name
                variable_field      = variable_fields{v};

                % This variable's true and fitted values
                variable_true       = True_Cylinder_Geometry.(variable_field);
                variable_mean       = Fitting_Method_Data.(variable_field).mu;
                variable_STD        = Fitting_Method_Data.(variable_field).sigma;
                variable_rel_STD    = Fitting_Method_Data.(variable_field).relative_sigma;
                variable_rel_error  = Fitting_Method_Data.(variable_field).relative_error;

                % Appended to the cells
                [variable_true_cell{v}, variable_mean_cell{v}, variable_rel_error_cell{v}, variable_STD_cell{v}, variable_rel_STD_cell{v}] = deal(variable_true, variable_mean, variable_rel_error, variable_STD, variable_rel_STD);
            end

            % The cells are converted to lists and then the table matrix
            [variable_true_list, variable_mean_list, variable_rel_error_list, variable_STD_list, variable_rel_STD_list] = deal(horzcat(variable_true_cell{:}), horzcat(variable_mean_cell{:}), horzcat(variable_rel_error_cell{:}), horzcat(variable_STD_cell{:}), horzcat(variable_rel_STD_cell{:}));
            table_data = [variable_true_list; variable_mean_list; variable_rel_error_list; variable_STD_list; variable_rel_STD_list]';

            % Data that is tabulated for each variable
            column_names = {'True [m]', sprintf('Fitted (%s) [m]', fitting_method), 'Rel. error [%]', 'Uncertainty [m]', 'Rel. uncert. [%]'};

            % The resulting table
            num_decimal_digits  = 3;
            output_format       = 'Float';
            table_file_name     = sprintf('Cylinder_Geometry_%s_Results.xls', fitting_label);
            results_table       = Table_Formatter(table_data, num_decimal_digits, output_format, variable_labels, column_names, table_file_name, Print);
                
            % The table is added to the structure and moved to the folder
            Results_Tables.(fitting_label) = results_table;
            
            movefile(table_file_name, save_folder_name);
        end

    %% Pairwise Student's t-test %%
        % The pairwise test is performed for each Monte Carlo sample of the fuzzy and least-squares geometry
        P_Values = Pairwise_Student_T_Test(MC_Geometry_Data.LS, MC_Geometry_Data.Optimal);
    
        P_values_cell = cell(1, number_variables);

        for v = 1 : number_variables
            variable_field      = variable_fields{v};
            P_values_list       = P_Values.(variable_field);
            P_values_cell{v}    = P_values_list;
        end

        P_values_list   = horzcat(P_values_cell{:});

        table_file_name = 'Cylinder_Geometry_T_Test.xls';
        row_name        = {'P [-] (t-test)'};
        output_format   = 'Exponential';
        t_test_table    = Table_Formatter(P_values_list, num_decimal_digits, output_format, row_name, variable_labels, table_file_name, Print);
        
        Results_Tables.t_test = t_test_table;

        movefile(table_file_name, save_folder_name);

    %% Plots %%  
        if Plot == true            
        %% Geometry overview with (weighted) point cloud %%
            %--% Preliminary %--%
            % Number of coordinates comprising the true and fitted cylinders
            number_coord    = 1e3;      

            % Colours for each scanner
            scanner_colours = cbrewer('qual', 'Set2', max(number_scanners, 3));     
            scanner_colours = max(scanner_colours, 0);
            scanner_colours = min(scanner_colours, 1);
            
            % Colour map for the point weights
            number_colours  = 1e3;
            weight_cmap     = cbrewer('seq', 'Greens', number_colours);
            weight_cmap     = max(weight_cmap, 0);
            weight_cmap     = min(weight_cmap, 1);

            [weight_min, weight_max]    = deal(0, max(point_weights_list));        % 0 is taken as the lower limit
            cmap_ind                    = round((number_colours - 1) * (point_weights_list - weight_min)/(weight_max - weight_min)) + 1;

            %--% Figure %--%
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            %----% Geometry on the first set of axes %----%  
            ax1 = axes;
            hold on
            grid on
            
            % The true cylinder surface
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(Cylinder_radius, Cylinder_length, Cylinder_centre, Cylinder_direction, number_coord);
            true_cyl_surf_name  = 'True cylinder';
            true_cyl_surf       = surf(ax1, cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', true_cyl_surf_name);

            % The true cylinder's caps
            pl_top = plot3(ax1, cylinder_coord_x(1, :), cylinder_coord_y(1, :), cylinder_coord_z(1, :), 'LineWidth', 2, 'color', 'k');
            pl_bot = plot3(ax1, cylinder_coord_x(2, :), cylinder_coord_y(2, :), cylinder_coord_z(2, :), 'LineWidth', 2, 'color', 'k');            
            pl_top.HandleVisibility = 'Off';
            pl_bot.HandleVisibility = 'Off';
        
            % The fitted cylinder surface
            Opt_GMM_radius      = Optim_Cylinder_Geometry_GMM.Cylinder_radius.mu;
            Opt_GMM_length      = Optim_Cylinder_Geometry_GMM.Cylinder_length.mu;
            Opt_GMM_centre      = Optim_Cylinder_Geometry_GMM.Cylinder_centre.mu;
            Opt_GMM_direction   = Optim_Cylinder_Geometry_GMM.Cylinder_direction.mu;

            [fitted_cylinder_coord_x, fitted_cylinder_coord_y, fitted_cylinder_coord_z, ~] = Cylinder_Surface_Generator(Opt_GMM_radius, Opt_GMM_length, Opt_GMM_centre, Opt_GMM_direction, number_coord);
            fit_cyl_surf_name   = 'Fitted cylinder';
            fit_cyl_surf        = surf(ax1, fitted_cylinder_coord_x, fitted_cylinder_coord_y, fitted_cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.10, 'LineWidth', 2, 'DisplayName', fit_cyl_surf_name);

            % The fitted cylinder's caps
            pl_top = plot3(ax1, fitted_cylinder_coord_x(1, :), fitted_cylinder_coord_y(1, :), fitted_cylinder_coord_z(1, :), 'LineWidth', 2, 'color', 'b');
            pl_bot = plot3(ax1, fitted_cylinder_coord_x(2, :), fitted_cylinder_coord_y(2, :), fitted_cylinder_coord_z(2, :), 'LineWidth', 2, 'color', 'b');            
            pl_top.HandleVisibility = 'Off';
            pl_bot.HandleVisibility = 'Off';
            
            % The scanner vectors
            pl_scan_vectors_name_cell   = cell(1, number_scanners);
            pl_scan_vectors_cell        = cell(1, number_scanners);

            for s = 1 : number_scanners
                % Distance from scanner to cylinder
                scanner_loc     = Scanner_loc_cell{s};
                scanner_vector  = scanner_loc - Cylinder_centre;
                range           = sqrt(sum(scanner_vector.^2));
                
                % The vector
                scanner_vector_scaled   = Cylinder_length * scanner_vector / norm(scanner_vector);
                scanner_string          = sprintf('Scanner %g, R = %.3g m', s, range);
                scanner_colour          = scanner_colours(s, :);

                pl_scan_vector = plot3(ax1, Cylinder_centre(1) + [0, scanner_vector_scaled(1)], Cylinder_centre(2) + [0, scanner_vector_scaled(2)], Cylinder_centre(3) + [0, scanner_vector_scaled(3)], 'LineWidth', 2, 'color', scanner_colour, 'DisplayName', scanner_string);

                % Added to the cell array for the legend
                pl_scan_vectors_name_cell{s}    = scanner_string;
                pl_scan_vectors_cell{s}         = pl_scan_vector;
            end

            %----% Point cloud shaded by weights on second set of axes %----%
            ax2 = axes;
            hold on
            grid on

            point_cloud_matrix  = vertcat(point_cloud_cell{:});
            point_cloud_sc_name = 'Point cloud';

            if strcmp(Point_Weighting, 'None')
                point_cloud_sc = scatter3(ax2, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'DisplayName', point_cloud_sc_name);
        
            else
                point_cloud_sc = scatter3(ax2, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 200, weight_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', point_cloud_sc_name);
        
                % Colorbar (on second axes)
                colormap(ax2, weight_cmap);
                cb = colorbar(ax2);
                shading interp
                clim([weight_min, weight_max])
                ylabel(cb, sprintf('%s [-]', Point_Weighting));
                cb.FontSize = 25;        
            end
    
            % Axes
            xlabel(ax1, 'x [m]')
            ylabel(ax1, 'y [m]')
            zlabel(ax1, 'z [m]')
    
            x_lim = [min(point_cloud_matrix(:, 1)) - Cylinder_radius, max(point_cloud_matrix(:, 1)) + Cylinder_radius];
            y_lim = [min(point_cloud_matrix(:, 2)) - Cylinder_radius, max(point_cloud_matrix(:, 2)) + Cylinder_radius];
            z_lim = [min(point_cloud_matrix(:, 3)) - Cylinder_radius, max(point_cloud_matrix(:, 3)) + Cylinder_radius];
    
            xlim([ax1, ax2], x_lim);
            ylim([ax1, ax2], y_lim);
            zlim([ax1, ax2], z_lim);
    
            ax1.DataAspectRatio = [1, 1, 1];
            ax2.DataAspectRatio = [1, 1, 1];
            view([ax1, ax2], 45, 45)

            % Legend
            legend(ax1, [true_cyl_surf, fit_cyl_surf, pl_scan_vectors_cell{:}, point_cloud_sc], [{true_cyl_surf_name}, {fit_cyl_surf_name}, pl_scan_vectors_name_cell(:)', {point_cloud_sc_name}], 'location', 'northoutside');
            
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
            
            % The figure is saved as .png and .fig and moved to the folder
            figure_name = 'Geometry_Overview';
            export_fig(1, [figure_name, '.fig']);
            export_fig(1, [figure_name, '.png']);
            
            movefile([figure_name, '*'], save_folder_name);
        
        %% Geometry parameter data plots %%              
            % Colours for each variable type            
            num_data_dim        = length(variable_labels);
            variable_colours    = cbrewer('qual', 'Set1', num_data_dim);
            variable_colours    = max(variable_colours, 0);
            variable_colours    = min(variable_colours, 1);

            % Data
            mean_GMM_field_selector     = @(fieldname) Optim_Cylinder_Geometry_GMM.(fieldname).mu;
            GMM_optim_geometry_cell     = cellfun(mean_GMM_field_selector, variable_fields, 'UniformOutput', false);
            GMM_optim_geometry          = horzcat(GMM_optim_geometry_cell{:});

            true_field_selector         = @(fieldname) True_Cylinder_Geometry.(fieldname);
            true_cylinder_geometry_cell = cellfun(true_field_selector, variable_fields, 'UniformOutput', false);
            true_cylinder_geometry      = horzcat(true_cylinder_geometry_cell{:});

            [MC_optimal_geometry_matrix, ~] = Structure_Data_Concatenation(MC_Geometry_Data.Optimal);
            [MC_LS_geometry_matrix, ~]      = Structure_Data_Concatenation(MC_Geometry_Data.LS);

            for p = 1 : num_data_dim
                % This data type's information
                variable_colour     = variable_colours(p, :);
                variable_type       = variable_labels{p};
                
                MC_opt_data     = MC_optimal_geometry_matrix(:, p);
                MC_LS_data      = MC_LS_geometry_matrix(:, p);
                true_value      = true_cylinder_geometry(p);
                fitted_value    = GMM_optim_geometry(p);
                                
                % This data type's figure
                figure(p + 10)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])   
                
                %--% Sorted data %--%
                subplot(1, 2, 1)
                hold on
                grid on
                
                sorted_MC_opt_data  = sort(MC_opt_data);
                sorted_MC_LS_data   = sort(MC_LS_data);
                plot(sorted_MC_opt_data, 'LineWidth', 2, 'color', variable_colour, 'DisplayName', sprintf('Optim. %s', variable_type));
                plot(sorted_MC_LS_data, 'LineWidth', 2, 'color', 'b', 'DisplayName', sprintf('LS %s', variable_type));
                plot([1, number_MC_iterations], [true_value, true_value], 'LineWidth', 2, 'LineStyle', '--', 'color', 'k', 'DisplayName', sprintf('%s true value', variable_type));
                plot([1, number_MC_iterations], [fitted_value, fitted_value], 'LineWidth', 2, 'LineStyle', ':', 'color', 'k', 'DisplayName', sprintf('%s fitted value', variable_type));
            
                % Axes
                ylabel(sprintf('%s [m]', variable_type));
                xlabel('Sorted sample [-]');
                xlim([1, max(number_MC_iterations, 2)]);         % In case the number of samples is 1

                % Legend
                legend('show', 'location', 'northoutside');
                
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off   
                
                %--% Histogram %--%
                subplot(1, 2, 2)
                hold on
                grid on
                
                % Histogram
                [number_bins, ~] = Histogram_Bins(MC_opt_data);

                hist_opt    = histogram(MC_opt_data, number_bins, 'FaceAlpha', 0.5, 'FaceColor', variable_colour, 'normalization', 'pdf', 'DisplayName', 'Monte Carlo optimiser samples');
                hist_LS     = histogram(MC_LS_data, number_bins, 'FaceAlpha', 0.5, 'FaceColor', 'b', 'normalization', 'pdf', 'DisplayName', 'Monte Carlo LS samples');
                f_max       = max([hist_opt.Values, hist_LS.Values]);       % The maximum density value

                % Fitted and true values
                plot([true_value, true_value], [0, f_max], 'LineWidth', 2, 'LineStyle', '--', 'color', 'k', 'DisplayName', sprintf('%s true value', variable_type));
                plot([fitted_value, fitted_value], [0, f_max], 'LineWidth', 2, 'LineStyle', ':', 'color', 'k', 'DisplayName', sprintf('%s fitted value', variable_type));

                % Axes
                xlabel([variable_type, ' [m]']);              
                ylabel('Probability density [1/m]');
                ylim([0, f_max]);
                
                % Legend
                legend('show', 'location', 'northoutside');

                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off    

                % The figure is saved as .png and .fig and moved to the folder
                figure_name = sprintf('Cylinder_Geometry_Data_%s', variable_type);
                export_fig(p + 10, [figure_name, '.fig']);
                export_fig(p + 10, [figure_name, '.png']);
            
                movefile([figure_name, '*'], save_folder_name);
            end

        %% General optimiser information %%
            % Mean objective values and its standard deviation for each step
            objective_steps_mean        = Objective_Steps.mu;
            objective_steps_STD         = Objective_Steps.sigma;
            objective_steps_STD_patch   = [objective_steps_mean + objective_steps_STD; flipud(objective_steps_mean - objective_steps_STD)];

            % Mean step sizes and its standard deviation for each step
            step_size_mean              = Step_Sizes.mu;
            step_size_STD               = Step_Sizes.sigma;
            step_size_STD_patch         = [step_size_mean + step_size_STD; flipud(step_size_mean - step_size_STD)];

            % Patch for the standard deviation
            max_num_optim_steps = max(number_steps_list);
            step_ind_patch      = [1 : max_num_optim_steps, max_num_optim_steps : -1 : 1]';

            figure(100)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   
            
            %--% Mean and standard deviation of objective value for each step %--%
            subplot(2, 3, 1)
            hold on
            grid on

            plot(1 : max_num_optim_steps, objective_steps_mean, 'LineWidth', 2, 'color', 'b', 'DisplayName', '\mu');
            patch(step_ind_patch, objective_steps_STD_patch, 'b', 'facealpha', 0.25, 'facecolor', 'b', 'edgecolor', 'none', 'DisplayName', '1 \sigma');
    
            xlabel('Number of steps [-]');
            ylabel('Q [-]');

            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Mean and standard deviation of the step size %--%
            subplot(2, 3, 2)
            hold on
            grid on

            plot(1 : max_num_optim_steps, step_size_mean, 'LineWidth', 2, 'color', 'c', 'DisplayName', '\mu');
            patch(step_ind_patch, step_size_STD_patch, 'c', 'facealpha', 0.25, 'facecolor', 'c', 'edgecolor', 'none', 'DisplayName', '1 \sigma');
    
            xlabel('Number of steps [-]');
            ylabel('\Delta [-]');

            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Objective value vs. number of steps %--%
            [objective_value_list_order, order] = sort(objective_value_list, 'ascend');
            number_steps_list_order             = number_steps_list(order);

            subplot(2, 3, 3)
            hold on
            grid on

            scatter(number_steps_list_order, objective_value_list_order, 'filled', 'MarkerFaceColor', 'b');

            xlabel('Number of steps [-]');
            ylabel('Q [-]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Histogram of number of steps %--%
            subplot(2, 3, 4)
            hold on
            grid on

            % Histogram
            max_number_steps    = max(number_steps_list);
            bin_edges           = 0.5 : 1 : max_number_steps + 0.5;

            histogram(number_steps_list, bin_edges, 'normalization', 'pdf');

            xlabel('Number of steps [-]');
            ylabel('f [-]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            %--% Histogram of the number of uncertainty updates %--%
            subplot(2, 3, 5)
            hold on
            grid on

            % Histogram
            [number_bins, ~] = Histogram_Bins(number_UU_iter_list);

            histogram(number_UU_iter_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', 'r', 'normalization', 'pdf');

            xlabel('Number of uncertainty updates [-]');
            ylabel('f [-]');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            % The figure is saved as .png and .fig and moved to the folder
            figure_name = 'Optimiser_Information';
            export_fig(100, [figure_name, '.fig']);
            export_fig(100, [figure_name, '.png']);
        
            movefile([figure_name, '*'], save_folder_name);            

        %% Plots for the optimisation geometry variables %%
            % Labels and units for the geometry variables used by the optimiser
            optim_geometry_labels       = Monte_Carlo_Loop.Optimiser_Diagnostics(1).parameter_labels;
            optim_geometry_units        = Monte_Carlo_Loop.Optimiser_Diagnostics(1).parameter_units;
            num_optim_geo_parameters    = length(optim_geometry_units);
    
            % The first and second derivatives of the objective w.r.t. the parameters at optimum
            First_Derivatives_cell          = {Monte_Carlo_Loop.Optimiser_Diagnostics.Optimum_First_Derivatives};
            First_Derivatives               = vertcat(First_Derivatives_cell{:});
            Norm_First_Derivatives_cell     = {Monte_Carlo_Loop.Optimiser_Diagnostics.Norm_Opt_First_Derivatives};
            Norm_First_Derivatives          = vertcat(Norm_First_Derivatives_cell{:});

            Second_Derivatives_cell         = {Monte_Carlo_Loop.Optimiser_Diagnostics.Optimum_Second_Derivatives};
            Second_Derivatives              = vertcat(Second_Derivatives_cell{:});
            Norm_Second_Derivatives_cell    = {Monte_Carlo_Loop.Optimiser_Diagnostics.Norm_Opt_Second_Derivatives};
            Norm_Second_Derivatives         = vertcat(Norm_Second_Derivatives_cell{:});

            % Colour map
            optim_colours    = cbrewer('qual', 'Set1', num_optim_geo_parameters);
            optim_colours    = max(optim_colours, 0);
            optim_colours    = min(optim_colours, 1);
    
            for g = 1 : num_optim_geo_parameters
                %--% This geometry parameter's data %--%
                geometry_label          = optim_geometry_labels{g};
                geometry_string         = strrep(geometry_label, '_', ' ');     % To prevent subscripts

                geometry_unit           = optim_geometry_units{g};
                geometry_colour         = optim_colours(g, :);
    
                geometry_steps_mean     = Geometry_Steps.(geometry_label).mu;
                geometry_steps_STD      = Geometry_Steps.(geometry_label).sigma;
                geometry_STD_patch      = [geometry_steps_mean + geometry_steps_STD; flipud(geometry_steps_mean - geometry_steps_STD)];
    
                geom_grad_steps_mean    = Geometry_Gradient_Steps.(geometry_label).mu;
                geom_grad_steps_STD     = Geometry_Gradient_Steps.(geometry_label).sigma;
                geom_grad_STD_patch     = [geom_grad_steps_mean + geom_grad_steps_STD; flipud(geom_grad_steps_mean - geom_grad_steps_STD)];
    
                geom_first_derivative_cell          = {First_Derivatives.(geometry_label)};
                geom_first_derivative_list          = vertcat(geom_first_derivative_cell{:});
                norm_geom_first_derivative_cell     = {Norm_First_Derivatives.(geometry_label)};
                norm_geom_first_derivative_list     = vertcat(norm_geom_first_derivative_cell{:});

                geom_second_derivative_cell         = {Second_Derivatives.(geometry_label)};
                geom_second_derivative_list         = vertcat(geom_second_derivative_cell{:});
                norm_geom_second_derivative_cell    = {Norm_Second_Derivatives.(geometry_label)};
                norm_geom_second_derivative_list    = vertcat(norm_geom_second_derivative_cell{:});

                %--% The plot %--%
                figure(200 + g)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])  
    
                % Parameter
                subplot(2, 3, 1)
                hold on
                grid on
    
                plot(1 : max_num_optim_steps, geometry_steps_mean, 'LineWidth', 2, 'color', geometry_colour, 'DisplayName', '\mu');
                patch(step_ind_patch, geometry_STD_patch, geometry_colour, 'facealpha', 0.25, 'facecolor', geometry_colour, 'edgecolor', 'none', 'DisplayName', '1 \sigma');
    
                xlabel('Number of steps [-]');
                ylabel(sprintf('%s [%s]', geometry_string, geometry_unit));
    
                legend('show', 'location', 'northoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
    
                % Histogram of first derivatives at optimum
                subplot(2, 3, 2)
                hold on
                grid on
    
                [number_bins, ~] = Histogram_Bins(geom_first_derivative_list);
    
                histogram(geom_first_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', geometry_colour, 'normalization', 'pdf');
    
                xlabel(sprintf('dQ/d%s [1/%s]', geometry_string, geometry_unit));
                ylabel('f [-]');
        
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % Histogram of normalised first derivatives at optimum
                subplot(2, 3, 3)
                hold on
                grid on
    
                [number_bins, ~] = Histogram_Bins(norm_geom_first_derivative_list);
    
                histogram(norm_geom_first_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', geometry_colour, 'normalization', 'pdf');
    
                xlabel(sprintf('dQ/d%s_n [-]', geometry_string));
                ylabel('f [-]');
        
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % Gradient of objective w.r.t. parameter
                subplot(2, 3, 4)
                hold on
                grid on
    
                plot(1 : max_num_optim_steps, geom_grad_steps_mean, 'LineWidth', 2, 'color', geometry_colour, 'DisplayName', '\mu');
                patch(step_ind_patch, geom_grad_STD_patch, geometry_colour, 'facealpha', 0.25, 'facecolor', geometry_colour, 'edgecolor', 'none', 'DisplayName', '1 \sigma');
    
                xlabel('Number of steps [-]');
                ylabel(sprintf('dQ/d%s [1/%s]', geometry_string, geometry_unit));
    
                legend('show', 'location', 'northoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off

                % Histogram of second derivatives at optimum
                subplot(2, 3, 5)
                hold on
                grid on
    
                [number_bins, ~] = Histogram_Bins(geom_second_derivative_list);
    
                histogram(geom_second_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', geometry_colour, 'normalization', 'pdf');
    
                xlabel(sprintf('d^2Q/d%s^2 [1/%s^2]', geometry_string, geometry_unit));
                ylabel('f [-]');
        
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                % Histogram of normalised second derivatives at optimum
                subplot(2, 3, 6)
                hold on
                grid on
    
                [number_bins, ~] = Histogram_Bins(norm_geom_second_derivative_list);
    
                histogram(norm_geom_second_derivative_list, number_bins, 'FaceAlpha', 1.0, 'FaceColor', geometry_colour, 'normalization', 'pdf');
    
                xlabel(sprintf('d^2Q/d%s_n^2 [-]', geometry_string));
                ylabel('f [-]');
        
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                % The figure is saved as .png and .fig and moved to the folder
                figure_name = sprintf('Optimiser_Steps_%s', geometry_label);
                export_fig(200 + g, [figure_name, '.fig']);
                export_fig(200 + g, [figure_name, '.png']);
            
                movefile([figure_name, '*'], save_folder_name);
            end  

        %% Point cloud (sampling) diagnostics %%
            % The expected Mahalanobis distance, propagation error and cylinder distance of each point for each iteration are analysed

            % Uncertainty and incidence angle of each point
            [sigma_radial_cell, sigma_prop_cell, incidence_angle_cell, ~, ~]    = Cylindrical_Object_Uncertainty(Cylinder_centre, Cylinder_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
            [sigma_radial_list, sigma_prop_list, PC_incidence_angle_list]          = deal(vertcat(sigma_radial_cell{:}), vertcat(sigma_prop_cell{:}), vertcat(incidence_angle_cell{:}));

            % Distributions
            alpha                       = 1 - Confidence_interval/100;
            Distr_Diagnostics           = false;
            Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

            % Expected cylinder distance (delta) and propagation error (epsilon) for the true geometry
            Bias_Plot                                           = false;
            [expected_delta_list, expected_epsilon_list, ~, ~]  = Cylinder_Expected_Bias(Cylinder_radius, Cylinder_centre, Cylinder_direction, Point_Cloud_Distributions, Scanner_loc_cell, Bias_Plot);
            
            % Put in the same order as the incidence angles of the unsampled point cloud
            [PC_incidence_angle_list, order]    = sort(PC_incidence_angle_list, 'ascend');
            PC_incidence_angle_list             = rad2deg(PC_incidence_angle_list);             % Conversion to degrees

            %--% Plot %--%
            % Data that is shown versus the incidence angle
            point_cloud_data_cell   = {expected_Mahal_dist_matrix, propagation_error_matrix, cylinder_distance_matrix};
            expected_data_cell      = {[], expected_epsilon_list, expected_delta_list};
            data_labels_cell        = {sprintf('E[M^%i]', distance_moment), sprintf('%s', '\epsilon'), sprintf('%s', '\delta')};
            data_units_cell         = {'-', 'm', 'm'};
            number_data_sets        = length(data_units_cell);
            data_cmap               = cbrewer('qual', 'Set2', max(number_data_sets, 3));

            % Incidence angles are transformed from radians to degrees and put in ascending order
            incidence_angle_list                    = rad2deg(incidence_angle_matrix(:));
            [incidence_angle_list, incidence_order] = sort(incidence_angle_list, 'ascend');

            % Plot
            figure(300)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])   

            % Subplots for each data set
            for d = 1 : number_data_sets
                % The data set
                data_list           = point_cloud_data_cell{d}(incidence_order);           % Sorted by the incidence angle 
                expected_data_list  = expected_data_cell{d};
                data_label          = data_labels_cell{d};
                data_unit           = data_units_cell{d};
                data_colour         = data_cmap(d, :);

                %--% Subplot %--%
                subplot(1, number_data_sets, d)
                hold on
                grid on

                % Data obtained during the Monte Carlo loop
                scatter(incidence_angle_list, data_list, 'MarkerFaceColor', data_colour, 'MarkerEdgeColor', 'none', 'DisplayName', 'Data');

                % Expected value (if available)
                if ~isempty(expected_data_list)
                    % Complete list
                    expected_data_list = expected_data_list(order);                                                                         % Put in the order of the incidence angles
                    plot(PC_incidence_angle_list, expected_data_list, 'LineWidth', 2, 'color', 'k', 'DisplayName', 'Expected value');       % Note that these are the incidence angles of the unsampled point cloud

                    % Mean value
                    mean_expected_value = mean(expected_data_list);
                    plot([0, max(PC_incidence_angle_list)], mean_expected_value*[1, 1], 'LineWidth', 2, 'color', 'k', 'LineStyle', ':', 'DisplayName', 'Mean expected value');
                end

                % Axes
                xlabel(sprintf('%s [deg]', '\alpha'));
                xlim([0, 90]);
                ylabel(sprintf('%s [%s]', data_label, data_unit));

                % Legend
                legend('show', 'location', 'northoutside');

                % Text
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);

                hold off
            end

            % The figure is saved as .png and .fig and moved to the folder
            figure_name = 'Sampled_Point_Cloud_Data';
            export_fig(300, [figure_name, '.fig']);
            export_fig(300, [figure_name, '.png']);
            
            movefile([figure_name, '*'], save_folder_name);

        %% Relationship between weight and geometry and uncertainty %%
            % The point weights are plotted w.r.t. the incidence angle, distance along the cylinder axis and uncertainty, both for the true and fitted geometry
            if ~strcmp(Point_Weighting, 'None')
                %--% Data %--%
                % The average fitted cylinder centre and direction (not taken from the GMM)
                fit_cyl_centre  = Fuzzy_MC_Geometry.Cylinder_centre.mu;
                fit_cyl_dir     = Fuzzy_MC_Geometry.Cylinder_direction.mu;

                % Data for the true geometry
                covariance_det_list     = 2*sigma_radial_list.^2 + sigma_prop_list.^2;

                [~, delta_list, ~]      = Point_to_Vector_Projection(point_cloud_matrix, Cylinder_direction, Cylinder_centre);
                cylinder_axis_dist_list = abs(delta_list);

                % Data for the fitted geometry
                [fit_sigma_radial_cell, fit_sigma_prop_cell, fit_incidence_angle_cell, ~, ~]    = Cylindrical_Object_Uncertainty(fit_cyl_centre, fit_cyl_dir, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
                [fit_sigma_radial_list, fit_sigma_prop_list, fit_incidence_angle_list]          = deal(vertcat(fit_sigma_radial_cell{:}), vertcat(fit_sigma_prop_cell{:}), vertcat(fit_incidence_angle_cell{:}));
                fit_covariance_det_list                                                         = 2*fit_sigma_radial_list.^2 + fit_sigma_prop_list.^2;

                [~, delta_list, ~]          = Point_to_Vector_Projection(point_cloud_matrix, fit_cyl_dir, fit_cyl_centre);
                fit_cylinder_axis_dist_list = abs(delta_list);

                % Collected into structures
                Geometry_Data       = struct('True_geometry', struct('covariance_determinant', covariance_det_list, 'incidence_angle', rad2deg(incidence_angle_list), 'cylinder_axis_distance', cylinder_axis_dist_list), ...                   % Note conversion to degrees
                                             'Fitted_geometry', struct('covariance_determinant', fit_covariance_det_list, 'incidence_angle', rad2deg(fit_incidence_angle_list), 'cylinder_axis_distance', fit_cylinder_axis_dist_list));        % Note conversion to degrees
                data_set_names      = fieldnames(Geometry_Data);
                number_data_sets    = length(data_set_names);

                Variable_Info       = struct('labels', struct('covariance_determinant', '|\Sigma|', 'incidence_angle', '\alpha', 'cylinder_axis_distance', '\delta'), ...
                                              'units', struct('covariance_determinant', 'm^2', 'incidence_angle', 'deg', 'cylinder_axis_distance', 'm'));
                variable_fields     = fieldnames(Variable_Info.labels);
                number_variables    = length(variable_fields);

                % Colours for each variable
                variable_cmap = cbrewer('qual', 'Set1', number_variables);

                %--% Plot %--%
                figure(400)
                % Set the size and white background color
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])   
            
                % Create the subplots
                Tiled_Chart_Layout_Main = tiledlayout(number_data_sets, 1);
                Tiled_Chart_Layout      = gobjects(1, number_data_sets);
                Axes                    = gobjects([number_data_sets, number_variables]);

                for d = 1 : number_data_sets
                    data_set_name = data_set_names{d};

                    Tiled_Chart_Layout(d) = tiledlayout(Tiled_Chart_Layout_Main, 1, number_variables);
                    Tiled_Chart_Layout(d).Layout.Tile = d;
                    
                    for v = 1 : number_variables
                        Axes(d, v) = nexttile(Tiled_Chart_Layout(d));

                        % The variable's data for this data set
                        variable_field  = variable_fields{v};
                        data_list       = Geometry_Data.(data_set_name).(variable_field);
                        variable_label  = Variable_Info.labels.(variable_field);
                        variable_unit   = Variable_Info.units.(variable_field);
                        variable_colour = variable_cmap(v, :);

                        % The correlation coefficient
                        corr_coeffs = corrcoef(point_weights_list, data_list);
                        corr_coeff  = corr_coeffs(1, 2);                        % The diagonals are auto-correlation

                        % Scatter plot
                        scatter(point_weights_list, data_list, 'MarkerFaceColor', variable_colour, 'MarkerEdgeColor', 'none', 'DisplayName', sprintf('%s = %.3g', '\rho', corr_coeff));

                        hold on
                        grid on
                        xlabel(sprintf('%s [-]', Point_Weighting));
                        ylabel(sprintf('%s [%s]', variable_label, variable_unit));

                        legend('show', 'location', 'east');
                    end

                    % The name of the geometry data set is used as title above the rows
                    row_title       = strrep(data_set_name, '_', ' ');
                    title(Tiled_Chart_Layout(d), row_title, 'FontSize', 15);
                end

                % Formatting
                set(Axes, 'FontSize', 15);
                set(Axes, 'LineWidth', 2);

                % The figure is saved
                figure_name = 'Point_Weighting_Diagnostics';

                export_fig(400, [figure_name, '.fig']);
                export_fig(400, [figure_name, '.png']);
            
                movefile([figure_name, '*'], save_folder_name);
            end
        end
end