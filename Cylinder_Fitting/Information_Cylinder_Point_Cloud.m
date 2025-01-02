% The information held in a point of the point cloud used for cylinder fitting is approximated as the delta in geometry when it is removed

function [Information_list, Relative_Deviation, Cylinder_Geometry_LOO] = Information_Cylinder_Point_Cloud(Parallel_Pool, Scanner_Parameters, Scanning_Parameters, Statistical_Values, Fitting_Parameters, Point_Cloud_Coord, Output_Decisions)

    %% Structure inputs %%
        % Parallel pool
        Parallel_Loop           = Parallel_Pool.Parallel_Loop;
        idle_timeout            = Parallel_Pool.idle_timeout;

        % Point cloud
        point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;
        number_points_list      = Point_Cloud_Coord.number_points_list;

        % Scanner parameters
        range_bias              = Scanner_Parameters.range_bias;

        % Scanning parameters
        Scanner_loc_cell        = Scanning_Parameters.Scanner_loc_cell;

        % Fitting parameters
        Distance_Computation    = Fitting_Parameters.Distance_Computation;

        % Statistical values
        Confidence_interval     = Statistical_Values.Confidence_interval;

        % Outputs
        Plot                    = Output_Decisions.Plot;
        Progress_Update         = Output_Decisions.Progress_Update;

    %% Manual inputs %%
        % Diagnostics
        Distr_Diagnostics       = false;        % [true, false] Shows how the distributions are created
        Mahal_Diagnostics       = false;        % [true, false] Diagnostics regarding E[M^2] computation

    %% Initiate the parallel pool
        if Parallel_Loop == true
            % The parallel pool is started
            number_cores_parallel = feature('numcores');                  % Checks the number of available cores
            Parallel_Pool_Starter(idle_timeout, number_cores_parallel);
        else
            % Otherwise the number of cores used in parallel are set to 0 to make the parfor run as a for
            number_cores_parallel = 0;
        end

    %% Leave-one-out cylinder fitting %%
        % Uniform weights are appended to the fitting parameters
        number_points                           = sum(number_points_list);
        Uniform_value_list                      = ones(number_points - 1, 1);      % Minus one as one point is always left out
        Fitting_Parameters.point_weights_list   = Uniform_value_list;

        % For the information estimate, the infinite cylinder geometry is saved
        Cylinder_Geometry_LOO               = struct('Cylinder_centre', [], 'Cylinder_direction', [], 'Cylinder_radius', []);      
        geometry_fields                     = fieldnames(Cylinder_Geometry_LOO);
        number_geometry_parameters          = length(geometry_fields);

        % Outputs within the loop are not desired
        Output_Decisions_LOO                = struct('Print', false, 'Plot', false, 'Diagnostics', false);

        % The cumulative number of beams is required to translate the point number to its scanner and beam indices
        cumul_beams_list                    = [0, cumsum(number_points_list)];
     
        % Progress counter
        DQ      = parallel.pool.DataQueue;
        tick    = 0;
        N       = number_points;
        afterEach(DQ, @ProgressUpdate);
    
        parfor (p = 1 : number_points, number_cores_parallel)
            % The associated scanner's beam
            scanner = find(p > cumul_beams_list, 1, 'last');
            beam    = p - cumul_beams_list(scanner);
    
            % It is removed from the point cloud and a new structure is created
            point_cloud_cell_LOO                    = point_cloud_cell;
            point_cloud_cell_LOO{scanner}(beam, :)  = [];

            number_points_list_LOO              = number_points_list;
            number_points_list_LOO(scanner)     = number_points_list_LOO(scanner) - 1;
    
            Point_Cloud_Coord_LOO = struct('point_cloud_cell', {point_cloud_cell_LOO}, 'number_points_list', number_points_list_LOO);
    
            % Cylinder fitting
            [Cylinder_Geometry_p, ~, ~, ~, ~] = Cylinder_Cross_Section_and_Direction_Estimation_fmincon(Point_Cloud_Coord_LOO, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Fitting_Parameters, Output_Decisions_LOO);
            
            % For consistency, the results are saved with positive sign in the third dimension
            Cylinder_direction_p                    = Cylinder_Geometry_p.Cylinder_direction;
            Cylinder_Geometry_p.Cylinder_direction  = sign(Cylinder_direction_p(3)) * Cylinder_direction_p;
            Cylinder_Geometry_LOO(p)                = Cylinder_Geometry_p;
    
            % Progress update
            if Progress_Update == true
                send(DQ, p);
            end
        end
            
    %% Geometry deviation %%
        % The relative deviation for each point
        Relative_Deviation      = struct();
        Avg_Cylinder_Geometry   = struct();
    
        % To remove scale-dependence and make the deviation relative, the values are normalised by the variables named here. If the field is empty, no division takes place
        normalisation_fields    = {'Cylinder_radius', '', 'Cylinder_radius'};          
    
        for g = 1 : number_geometry_parameters
            % This parameter's data and its mean
            geometry_parameter  = geometry_fields{g};            
            geometry_cell       = {Cylinder_Geometry_LOO(:).(geometry_parameter)};
            geometry_matrix     = vertcat(geometry_cell{:});
    
            mean_geometry                               = mean(geometry_matrix, 1);
            Avg_Cylinder_Geometry.(geometry_parameter)  = mean_geometry;
    
            % The one by which it is normalised
            normalis_parameter  = normalisation_fields{g};
    
            if ~isempty(normalis_parameter)
                normalising_cell    = {Cylinder_Geometry_LOO(:).(normalis_parameter)};
                normalising_matrix  = vertcat(normalising_cell{:});
                normalising_value   = mean(normalising_matrix, 1);
            else
                normalising_value   = 1;
            end
    
            % The relative deviation
            relative_dev_matrix                     = (mean_geometry - geometry_matrix) / normalising_value * 100;      % Note the conversion to percent
            Relative_Deviation.(geometry_parameter) = relative_dev_matrix;
        end

        % The total relative deviation for each left out point is the information metric
        total_relative_dev_cell     = struct2cell(Relative_Deviation);
        total_relative_dev_matrix   = horzcat(total_relative_dev_cell{:});
        total_relative_dev_list     = sqrt(sum(total_relative_dev_matrix.^2, 2));

    %% Expected Mahalanobis distance %%
        % As anomalous points are also expected to result in a large deviation, the points are weighted by the expected Mahalanobis distance as well
        % The average geometry estimate during LOO computation is used
        avg_cyl_dir         = Avg_Cylinder_Geometry.Cylinder_direction;
        avg_cyl_radius      = Avg_Cylinder_Geometry.Cylinder_radius;
        avg_cyl_centre      = Avg_Cylinder_Geometry.Cylinder_centre;

        % The uncertainty of each point, given the geometry
        [sigma_radial_cell, sigma_prop_cell, ~, ~] = Cylindrical_Object_Uncertainty(avg_cyl_centre, avg_cyl_dir, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);

        alpha                       = 1 - Confidence_interval/100;
        Point_Cloud_Distributions   = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Distr_Diagnostics);

        % The expected squared Mahalanobis distance
        distance_moment = 2;

        if strcmp(Distance_Computation, 'Line_Approx')
            [~, expected_Mahal_distance_list] = Expected_Mahal_Distance_Cylinder_Line_Approx(avg_cyl_centre, avg_cyl_radius, avg_cyl_dir, distance_moment, Point_Cloud_Distributions, Scanner_loc_cell, Mahal_Diagnostics);
        elseif strcmp(Distance_Computation, 'Numerical')
            [~, expected_Mahal_distance_list] = Expected_Mahalanobis_Distance_Cylinder_Numerical(avg_cyl_centre, avg_cyl_radius, avg_cyl_dir, distance_moment, Point_Cloud_Distributions, Statistical_Values, Mahal_Diagnostics);
        end

    %% Information %%
        % The information is taken to be the total relative deviation divided by the square root of the expected squared Mahalanobis distance
        Information_list = total_relative_dev_list ./ sqrt(expected_Mahal_distance_list);

        % The values are normalised s.t. the average is 1
        Information_list = Information_list / mean(Information_list);

    %% Plots %%
        if Plot == true
            % All the data sets that will be plotted
            num_geo_dim         = size(total_relative_dev_matrix, 2);
            relative_dev_cell   = mat2cell(total_relative_dev_matrix, number_points, ones(1, num_geo_dim));
            
            data_sets           = [{Information_list}, {expected_Mahal_distance_list}, {total_relative_dev_list}, relative_dev_cell];
            num_data_sets       = length(data_sets);
    
            % Their labels and colour maps
            data_set_labels = {'I', 'E[M^2]', '\Delta G', '\Delta c_x', '\Delta c_y', '\Delta c_z', '\Delta v_x', '\Delta v_y', '\Delta v_z', '\Delta r'};
            data_set_units  = {'-', '-', '%', '%', '%', '%', '%', '%', '%', '%'};
            data_set_cmaps  = {'Blues', 'Reds', 'Greys', 'PiYG', 'PiYG', 'PiYG', 'RdBu', 'RdBu', 'RdBu', 'PRGn'};           % Each category gets its colour map
            cmap_categories = {'seq', 'seq', 'seq', 'div', 'div', 'div', 'div', 'div', 'div', 'div'};                       % Sequential for purely positive ones, divergent for negative-positive
            num_colours     = 1e3;
        
            % The length is estimated from the point cloud
            point_cloud_matrix  = vertcat(point_cloud_cell{:});
            delta_point_matrix  = point_cloud_matrix - mean(point_cloud_matrix, 1);
            point_norms         = sqrt(sum(delta_point_matrix.^2, 2));
            cyl_length          = 2*max(point_norms);
    
            num_points_cyl                                              = 1e3;
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~]   = Cylinder_Surface_Generator(avg_cyl_radius, cyl_length, avg_cyl_centre, avg_cyl_dir, num_points_cyl);
                
            % The resulting plots
            link_struct = struct('link', []);       % The links need to be saved for linkprop to remain active
    
            for d = 1 : num_data_sets
                % This data set's data and label
                data_set    = data_sets{d};
                data_label  = data_set_labels{d};
                data_unit   = data_set_units{d};
    
                % The colour map
                cmap_name   = data_set_cmaps{d};
                cmap_cat    = cmap_categories{d};
    
                data_cmap   = cbrewer(cmap_cat, cmap_name, num_colours);
                data_cmap   = max(data_cmap, 0);
                data_cmap   = min(data_cmap, 1);
    
                if strcmp(cmap_cat, 'div')                  % To ensure that the diverging colour map is appropriately centered, the limits are equal in both directions
                    data_max_abs    = max(abs(data_set));
                    data_min        = -data_max_abs;
                    data_max        = data_max_abs;
                else
                    data_min        = min(data_set);
                    data_max        = max(data_set);
                end
    
                if data_min == 0 && data_max == 0           % In case there was zero change for any point
                    data_max = 1e-6;
                end
    
                cmap_ind = round((num_colours - 1) * (data_set - data_min)/(data_max - data_min)) + 1;  
    
                figure(d);
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
                pt_name = 'Point cloud';
                pt_sc   = scatter3(ax2, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 200, data_cmap(cmap_ind, :), 'Marker', '.', 'DisplayName', pt_name);
        
                % Colorbar (on second axes)
                colormap(ax2, data_cmap);
        
                cb = colorbar(ax2);
                shading interp
                clim([data_min, data_max])
                ylabel(cb, sprintf('%s [%s]', data_label, data_unit));
                cb.FontSize = 25;        
        
                % Axes
                xlabel(ax1, 'x [m]')
                ylabel(ax1, 'y [m]')
                zlabel(ax1, 'z [m]')
        
                x_lim = [min(cylinder_coord_x(:)) - avg_cyl_radius, max(cylinder_coord_x(:)) + avg_cyl_radius];
                y_lim = [min(cylinder_coord_y(:)) - avg_cyl_radius, max(cylinder_coord_y(:)) + avg_cyl_radius];
                z_lim = [min(cylinder_coord_z(:)) - avg_cyl_radius, max(cylinder_coord_z(:)) + avg_cyl_radius];
        
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
                axis_dimensions = [ax1.Position; ax2.Position];
                axis_starts     = max(axis_dimensions(:, 1:2), [], 1);
                axis_sizes      = min(axis_dimensions(:, 3:4), [], 1);    
                ax1.Position    = [axis_starts, axis_sizes];
                ax2.Position    = [axis_starts, axis_sizes];
                ax2.Visible     = 'off';
        
                link                = linkprop([ax1, ax2], {'View', 'XLim', 'YLim', 'ZLim'});
                link_struct(d).link = link;
        
                hold([ax1, ax2], 'off');        
    
                % The figure is saved
                figure_name = sprintf('Point_Cloud_%s.png', data_label);
                figure_name = strrep(figure_name, '\', '');
                figure_name = strrep(figure_name, ' ', '_');
               
                export_fig(d, figure_name);
            end
    
            % Pause message
            disp('The information plots have been created and saved. The script continues automatically.');
            close(1 : num_data_sets);
        end

    %% Information progress function %%
        function ProgressUpdate(~)
            % Ensures that at most every percent is printed
            tick = tick + 1;    

            progress_last   = floor((tick - 1) / N * 100);
            progress        = floor(tick / N * 100);

            if progress - progress_last >= 1
                fprintf('   Information assessment progress: %g %% \n', progress);
            end            
        end

end