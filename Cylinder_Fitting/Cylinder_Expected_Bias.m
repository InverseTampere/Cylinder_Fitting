% The expected values (and thus bias) in the cylinder distance delta and propagation error epsilon are computed here for the given point cloud and cylinder

function [expected_delta_list, expected_epsilon_list, expected_delta_cell, expected_epsilon_cell] = Cylinder_Expected_Bias(cylinder_radius, cylinder_centre, cylinder_direction, Point_Cloud_Distributions, Scanner_loc_cell, Plot)

    %% Structure inputs %%
        % Point cloud distributions
        distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
        number_distributions_list   = Point_Cloud_Distributions.number_distributions_list;
        total_num_distributions     = Point_Cloud_Distributions.number_distributions;

    %% Manual inputs %%
        % Sampling the distributions
    	number_samples 		    = 1e2;				% [-] Samples taken for each dimension
	    number_STD 			    = 03;				% [-] Number of standard deviations used to space the samples

        % Numerical margin
        numerical_margin        = 1e-2;             % [-] To detect whether or not the projected propagation vector is correct

        % Diagnostics
        Proj_Diagnostics        = false;            % [true, false] Projection of the distributions onto the cylinder's cross-section

    %% Point cloud projection onto cross-section %%
        % To simplify the equations, coordinate bases are defined individually for each scanner
        [~, ~, ~, projection_vector_base_cell] = Vector_Basis_Cylinder_Cross_Section_Projection(cylinder_centre, cylinder_direction, Scanner_loc_cell);

        % The propagation vector are retrieved from the distribution axes, which should always be the last
        Prop_Vector_Dim_fun = @(distr_axes) distr_axes(end, :);    
        prop_vector_cell    = cellfun(Prop_Vector_Dim_fun, distribution_axes_cell, 'UniformOutput', false);
        prop_vector_matrix  = vertcat(prop_vector_cell{:});
        num_dim             = size(prop_vector_matrix, 2);

        % The transformed properties are stored in matrix format
        cumul_beams_list    = [0, cumsum(number_distributions_list)];
        number_scanners     = length(Scanner_loc_cell);

        proj_mu_cell            = cell(total_num_distributions, 1);
        proj_sigmae_cell        = cell(total_num_distributions, 1);
        proj_prop_vector_cell   = cell(total_num_distributions, 1);

        scaling_factor_list     = zeros(total_num_distributions, 1);        % To scale the lengths back from 2D to 3D

        for s = 1 : number_scanners
            % This scanner's indices
            ind_first   = cumul_beams_list(s) + 1;
            ind_last    = cumul_beams_list(s + 1);

            % The vector basis
            projection_vector_basis = projection_vector_base_cell{s};       % Vector basis s.t. the scanner lies on the y-z plane and the cylinder axis is the z-axis

            % Distributions projected on the cross-sectional plane
            Projected_Distributions_scanner = Multivariate_Normal_Plane_Projection(projection_vector_basis, cylinder_centre, Point_Cloud_Distributions, Proj_Diagnostics);

            % Projected mu, sigmae and propagation vector in the cross-sectional coordinate frame
            mu_Q_cell           = Projected_Distributions_scanner.Plane.Projection.mu(ind_first : ind_last);
            mu_Q_matrix         = vertcat(mu_Q_cell{:});
            sigmae_Q_cell       = Projected_Distributions_scanner.Plane.Projection.sigmae(ind_first : ind_last);
            sigmae_Q_matrix     = vertcat(sigmae_Q_cell{:});
            distr_axes_Q_cell   = Projected_Distributions_scanner.Plane.Projection.distr_axes(ind_first : ind_last);

            % The propagation vectors projected onto the plane
            prop_vector_Q_3D_matrix     = (projection_vector_basis * prop_vector_matrix(ind_first : ind_last, :)')';
            prop_vector_Q_3D_matrix     = prop_vector_Q_3D_matrix(:, 1 : num_dim - 1);      % The vector component orthogonal to the plane is always zero and therefore removed
            prop_vector_3D_norm_list    = sqrt(sum(prop_vector_Q_3D_matrix.^2, 2));
            prop_vector_Q_3D_matrix     = prop_vector_Q_3D_matrix ./ prop_vector_3D_norm_list;

            scaling_factor_list(ind_first : ind_last) = prop_vector_3D_norm_list;

            % The propagation vector is ensured to be last by checking its dot product with the projected propagation vector
            number_distributions    = number_distributions_list(s);
            dot_product_matrix      = zeros(number_distributions, 2);    

            for d = 1 : 2
                % This dimension's axes
                Vector_Dim_fun      = @(vector_matrix) vector_matrix(d, :);
                
                distr_axes_d_cell   = cellfun(Vector_Dim_fun, distr_axes_Q_cell, 'UniformOutput', false);
                distr_axes_d_matrix = vertcat(distr_axes_d_cell{:});

                % The dot products with the propagation vectors
                dot_product_list            = dot(prop_vector_Q_3D_matrix, distr_axes_d_matrix, 2);
                dot_product_matrix(:, d)    = dot_product_list;
            end

            % Propagation vectors are those with an absolute dot product of 1
            prop_vector_bool        = abs(dot_product_matrix) > 1 - numerical_margin;
            prop_vector_bool_cell   = mat2cell(prop_vector_bool, ones(1, number_distributions), num_dim - 1);

            Matrix_Bool_fun         = @(row_boolean, matrix) matrix(row_boolean, :);
            prop_vector_Q_cell      = cellfun(Matrix_Bool_fun, prop_vector_bool_cell, distr_axes_Q_cell, 'UniformOutput', false);
            prop_vector_Q_matrix    = vertcat(prop_vector_Q_cell{:});

            % Check to ensure that at least one of the vectors is found to be the propagation vector
            if sum(prop_vector_bool, 'all') ~= number_distributions
                warning('The propagation vector was not accurately detected for all distributions');
            end

            % If the propagation vector is in the first column, the uncertainties need to be switched
            dim_switch_bool                     = prop_vector_bool(:, 1);
            sigmae_Q_matrix(dim_switch_bool, :) = fliplr(sigmae_Q_matrix(dim_switch_bool, :));

            % If the sign of the dot product is -1, that means the signs need to be reversed
            dot_product_matrix_T    = dot_product_matrix';
            prop_vector_dot_product = dot_product_matrix_T(prop_vector_bool');
            sign_list               = sign(prop_vector_dot_product);
            prop_vector_Q_matrix    = sign_list .* prop_vector_Q_matrix;

            % Converted to cell arrays and inserted
            mu_Q_cell           = mat2cell(mu_Q_matrix, ones(1, number_distributions), num_dim - 1);
            sigmae_Q_cell       = mat2cell(sigmae_Q_matrix, ones(1, number_distributions), num_dim - 1);
            prop_vector_Q_cell  = mat2cell(prop_vector_Q_matrix, ones(1, number_distributions), num_dim - 1);

            proj_mu_cell(ind_first : ind_last)          = mu_Q_cell;
            proj_sigmae_cell(ind_first : ind_last)      = sigmae_Q_cell;
            proj_prop_vector_cell(ind_first : ind_last) = prop_vector_Q_cell;
        end

    %% Expected cylinder distance %%
        % The expected cylinder distance is equal in 2D and 3D
        Expected_Delta_fun  = @(distr_sigmae, prop_vector, mu) Expected_Delta(distr_sigmae, prop_vector, mu, cylinder_radius, number_samples, number_STD);
        expected_delta_list = cellfun(Expected_Delta_fun, proj_sigmae_cell, proj_prop_vector_cell, proj_mu_cell);

        % Cell array for each scanner
        expected_delta_cell = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's indices
            ind_first   = cumul_beams_list(s) + 1;
            ind_last    = cumul_beams_list(s + 1);

            % Its expected delta values
            expected_delta_cell{s} = expected_delta_list(ind_first : ind_last);
        end

    %% Expected propagation error %%
        % The expected propagation error in 2D
        Expected_Epsilon_fun        = @(distr_sigmae, prop_vector, mu) Expected_Epsilon(distr_sigmae, prop_vector, mu, cylinder_radius, number_samples, number_STD);
        expected_epsilon_2D_list    = cellfun(Expected_Epsilon_fun, proj_sigmae_cell, proj_prop_vector_cell, proj_mu_cell);

        % The scaling factor to 3D comes from projecting the 2D propagation vectors onto the 3D one
        expected_epsilon_list = scaling_factor_list .* expected_epsilon_2D_list;

        % Cell array for each scanner
        expected_epsilon_cell = cell(1, number_scanners);

        for s = 1 : number_scanners
            % This scanner's indices
            ind_first   = cumul_beams_list(s) + 1;
            ind_last    = cumul_beams_list(s + 1);

            % Its expected epsilon values
            expected_epsilon_cell{s} = expected_epsilon_list(ind_first : ind_last);
        end

    %% Plots %%
        if Plot == true
            % The data sets to be plotted
            data_set_cell       = {expected_delta_list * 1e3, expected_epsilon_list * 1e3};     % Note conversion from m to mm
            data_set_variables  = {'E[\delta]', 'E[\epsilon]'};
            data_set_units      = {'mm', 'mm'};
            number_data_sets    = length(data_set_units);

            number_colours  = 1e3;
            data_cmap_cell  = Colorbrewer_Colour_Maps('Diverging', number_data_sets, number_colours);

            % The point cloud
            point_cloud_matrix  = vertcat(distribution_mu_cell{:});

            % Estimate of the cylinder length
            [~, delta_list, ~]  = Point_to_Vector_Projection(point_cloud_matrix, cylinder_direction, cylinder_centre);
            cylinder_length     = 2*max(abs(delta_list));

            for d = 1: number_data_sets
                % The data
                data_list       = data_set_cell{d};
                data_variable   = data_set_variables{d};
                data_unit       = data_set_units{d};
                data_cmap       = data_cmap_cell{d};
                
                % Colours for the data points
                data_lim            = max(abs(data_list));
                colour_ind          = round((number_colours - 1) * (data_list + data_lim)/(2*data_lim)) + 1;
                [colour_ind, order] = sort(colour_ind);
                data_colours        = data_cmap(colour_ind, :);

                %--% Plot %--%
                figure(d)
                % Size and white background
                set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
                set(gcf, 'color', [1, 1, 1])    

                for s = 1 : number_scanners
                    % This scanner's data
                    ind_first   = cumul_beams_list(s) + 1;
                    ind_last    = cumul_beams_list(s + 1);

                    point_cloud_matrix_s    = point_cloud_matrix(ind_first : ind_last, :);
                    scanner_loc             = Scanner_loc_cell{s};
                    data_colours_s          = data_colours(ind_first : ind_last, :);

                    % Estimate of the range
                    point_cloud_matrix_sc   = point_cloud_matrix_s - scanner_loc;
                    range_list              = sqrt(sum(point_cloud_matrix_sc.^2, 2));
                    scanner_range           = round(mean(range_list));

                    %--% Subplot %--%
                    % First set of axes for the cylinder
                    ax1 = axes;
                    hold on
                    grid on
            
                    number_coord = 1e3;
                    [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(cylinder_radius, cylinder_length, cylinder_centre, cylinder_direction, number_coord);
                    surf(ax1, cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.25, 'LineWidth', 2);

                    % Second set of axes for the point cloud coloured by the data set
                    ax2 = axes;
                    hold on
                    grid on
                    scatter3(ax2, point_cloud_matrix(order, 1), point_cloud_matrix(order, 2), point_cloud_matrix(order, 3), 1e2, data_colours_s, 'filled');

                    % The colourmap is applied to the second set of axes for the data
                    colormap(ax2, data_cmap);
                    cb = colorbar(ax2);
                    shading interp
                    clim(data_lim * [-1, 1])
            
                    colortitlehandle    = get(cb, 'Title');
                    titlestring         = sprintf('%s [%s]', data_variable, data_unit);
                    set(colortitlehandle, 'String', titlestring);
                    cb.FontSize         = 15;
            
                    % Text        
                    xlabel(ax1, 'x [m]');
                    ylabel(ax1, 'y [m]');
                    zlabel(ax1, 'z [m]');
            
                    set(ax1, 'FontSize', 15);
                    set(ax1, 'LineWidth', 2);  
                                    
                    % Axes dimensions      
                    x_lim = [min(point_cloud_matrix(:, 1)), max(point_cloud_matrix(:, 1))] + cylinder_radius*[-1, 1];
                    y_lim = [min(point_cloud_matrix(:, 2)), max(point_cloud_matrix(:, 2))] + cylinder_radius*[-1, 1];
                    z_lim = [min(point_cloud_matrix(:, 3)), max(point_cloud_matrix(:, 3))] + cylinder_radius*[-1, 1];
            
                    xlim([ax1, ax2], x_lim);
                    ylim([ax1, ax2], y_lim);
                    zlim([ax1, ax2], z_lim);
            
                    ax1.DataAspectRatio = [1, 1, 1];
                    ax2.DataAspectRatio = [1, 1, 1];
    
                    % Move to subplot
                    subplot(1, number_scanners, s, ax1);
                    subplot(1, number_scanners, s, ax2);
            
                    title_str = sprintf('Scanner %i, R = %i m', s, scanner_range);
                    title(ax1, title_str);
    
                    % Change the size, link the axes and make the second set invisible
                    axis_dimensions = [ax1.Position; ax2.Position];
                    axis_starts     = max(axis_dimensions(:, 1:2), [], 1);
                    axis_sizes      = min(axis_dimensions(:, 3:4), [], 1);    
                    ax1.Position    = [axis_starts, axis_sizes];
                    ax2.Position    = [axis_starts, axis_sizes];
                    ax2.Visible     = 'off';
            
                    view([ax1, ax2], -45, 45);
                    linkaxes([ax1, ax2]);
                    linkprop([ax2 ax1], {'View', 'XLim', 'YLim', 'ZLim'});
                end
            end

            % Pause message
            disp('The expected bias has been computed. The figure will close and script end upon a key-press.');
            pause();

            close(1 : number_data_sets);
        end

    %% Local functions %%
        % The inputs are expected to have been projected onto the cross-sectional plane, centered on the circle's centre (i.e. cylinder axis)
        % a and b are the radial and propagation axes of the distribution respectively

	    % Circle distance
        function E_delta = Expected_Delta(distr_sigmae, vector_b, mu, circle_radius, number_samples, number_STD)
            % The circle centre in the distribution's coordinate frame
            c_a = -vector_b(2)*mu(1) + vector_b(1)*mu(2);
            c_b = -vector_b(1)*mu(1) - vector_b(2)*mu(2);
    
		    % Full array of samples 
            [sigma_a, sigma_b]  = deal(distr_sigmae(1), distr_sigmae(2));
		    sample_a_list       = sigma_a*number_STD*linspace(-1, 1, number_samples);
		    sample_b_list       = sigma_b*number_STD*linspace(-1, 1, number_samples);
		    
		    [sample_a_matrix, sample_b_matrix] 	= meshgrid(sample_a_list, sample_b_list);
		    sample_a_list 						= sample_a_matrix(:);
		    sample_b_list 						= sample_b_matrix(:);        
    
            % Distance from each sample to the circle centre
            sample_circle_dist_list = sqrt((sample_a_list - c_a).^2 + (sample_b_list - c_b).^2);
    
            % Probability density of each sample
            prob_density_list       = 1/(2*pi*sigma_a*sigma_b) * exp(-1/2*sample_a_list.^2/sigma_a^2 - 1/2*sample_b_list.^2/sigma_b^2);
    
            % Resulting circle distance
            E_delta = -circle_radius + sum(sample_circle_dist_list .* prob_density_list) / sum(prob_density_list);
        end

        % Propagation error
        function E_epsilon = Expected_Epsilon(distr_sigmae, vector_b, mu, circle_radius, number_samples, number_STD)
            % The circle centre in the distribution's coordinate frame
            c_a = -vector_b(2)*mu(1) + vector_b(1)*mu(2);
            c_b = -vector_b(1)*mu(1) - vector_b(2)*mu(2);

		    % Full array of samples 
            [sigma_a, sigma_b]  = deal(distr_sigmae(1), distr_sigmae(2));
		    sample_a_list       = sigma_a*number_STD*linspace(-1, 1, number_samples);
		    sample_b_list       = sigma_b*number_STD*linspace(-1, 1, number_samples);

            % First intersections of each sample with the circle in the propagation direction
            points_b_list           = repmat(c_b - 2*circle_radius, [number_samples, 1]);       % Such that the first intersection is always on the correct side of the circle
            points_matrix           = [sample_a_list', points_b_list];
            vector_matrix           = repmat([0, 1], [number_samples, 1]);                      % The propagation vector has been defined as the y-axis
            circle_centre           = [c_a, c_b];
            Intersect_Diagnostics   = false;

            [first_intersects_matrix, ~, ~, intersection_bool]  = Vector_Circle_Intersection(points_matrix, vector_matrix, circle_centre, circle_radius, Intersect_Diagnostics);
            first_intersect_b_list                              = first_intersects_matrix(:, 2);

            % Mesh of samples which have an intersection
            sample_a_list                       = sample_a_list(intersection_bool);
		    [sample_a_matrix, sample_b_matrix] 	= meshgrid(sample_a_list, sample_b_list);

            % Propagation error of each sample
            prop_error_matrix   = sample_b_matrix - first_intersect_b_list';

            % Probability density of each sample which had an intersect
            prob_density_matrix = 1/(2*pi*sigma_a*sigma_b) * exp(-1/2*sample_a_matrix.^2/sigma_a^2 - 1/2*sample_b_matrix.^2/sigma_b^2);
    
            % Resulting expected propagation error
            E_epsilon = sum(prop_error_matrix .* prob_density_matrix, 'all') / sum(prob_density_matrix, 'all');

            % If there are no intersections, an expected value of 0 is given
            if sum(intersection_bool) == 0
                E_epsilon = 0;
            end
        end	
end