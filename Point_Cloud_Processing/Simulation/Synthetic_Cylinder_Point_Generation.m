% This script generates synthetic points based on the given laser beam uncertainty and scanning geometry, 
% and generates 3D multivariate normal distributions based on them.
% Note that the scanning vectors are spaced uniformly and parallel to the Cartesian axes.

% Printed statements giving a basic overview of the point cloud are shown if Print is true.

function [Point_Cloud_Coord, Point_Cloud_Distributions, Uncertainty_Coverage] = Synthetic_Cylinder_Point_Generation(True_Cylinder_Geometry, Scanning_Parameters, Scanner_Parameters, Statistical_Values, Output_Decisions)

    %% Inputs %%   
        % Scanned geometry
        Cylinder_centre     = True_Cylinder_Geometry.Cylinder_centre;
        Cylinder_direction  = True_Cylinder_Geometry.Cylinder_direction;
        Cylinder_length     = True_Cylinder_Geometry.Cylinder_length;
        Cylinder_radius     = True_Cylinder_Geometry.Cylinder_radius;
    
        % Scanning parameters
        Initial_Uncertainty = Scanning_Parameters.Initial_Uncertainty;
        Scanner_loc_cell    = Scanning_Parameters.Scanner_loc_cell;
        number_scanners     = Scanning_Parameters.number_scanners;
        
        % Scanner parameters
        laser_resolution    = Scanner_Parameters.laser_resolution;                 
        range_bias          = Scanner_Parameters.range_bias;                      
        max_incidence_angle = Scanner_Parameters.max_incidence_angle;       
        angular_accuracy    = Scanner_Parameters.angular_accuracy;   
        beam_divergence     = Scanner_Parameters.beam_divergence;

        % Statistical values
        Confidence_interval = Statistical_Values.Confidence_interval;
        
        % Output decisions
        Print               = Output_Decisions.Print;
        
    %% Manual inputs %%
        % Intermediate diagnostic plots may be shown as well if desired
        Diagnostics_Distr_Generation    = false;        % [true, false] Produces a plot of the multivariate normal distributions
        Diagnostics_Projection          = false;        % [true, false] Shows a plot of the projected distributions
        Sampling_Diagnostics            = false;        % [true, false] Sampled point cloud
        Diagnostics                     = false;        % [true, false] The final generated point cloud

        % The coverage angle is determined to within the specified accuracy
        coverage_resolution             = 0.1;          % [deg]
        
    %% The scanner beam vectors %%
        %--% The limiting z-coordinates of the top and bottom of the cylinder %--%
        % The normalised vector of the cylinder
        Cylinder_direction  = Cylinder_direction ./ sqrt(sum(Cylinder_direction.^2, 2));
        
        % Elevation angle
        num_dim             = length(Cylinder_centre);
        Cylinder_elev_angle = atan(Cylinder_direction(num_dim) / norm(Cylinder_direction(1 : num_dim - 1)));
        
        % Centres of the top and bottom 'caps' (including margin)
        cap_centres = Cylinder_centre + [1; -1]*(Cylinder_length/2*Cylinder_direction + Cylinder_radius*cos(Cylinder_elev_angle));        % Top; Bottom
        
        % The largest extent of a slice in the x-y plane
        if Cylinder_direction(1) == 1 || Cylinder_direction(2) == 1         % The cylinder axis is parallel to the plane
            Cyl_slice_extent    = Cylinder_length;
        else
            gamma               = pi/2 - atan(Cylinder_direction(3) / Cylinder_direction(1)); 
            Cyl_slice_extent    = 2*Cylinder_radius * sec(gamma);
        end
        
        %--% The scanning beam vectors %--% 
        % Note that the grid is parallel to the Cartesian axes, and thus independent of cylinder orientation
        beam_vector_cell = cell(1, number_scanners);
            
        for s =  1 : number_scanners
            % Scanner location
            Scanner_loc = Scanner_loc_cell{s};

            % Vector from the scanner to the cylinder centre on the x-y-plane
            centre_vector           = Cylinder_centre - Scanner_loc;
            centre_vector(num_dim)  = 0;
            centre_vector           = centre_vector / norm(centre_vector);

            % Rotation matrix s.t. this vector becomes the x-axis
            rotation_matrix = [centre_vector; -centre_vector(2), centre_vector(1), 0; 0, 0, 1];

            % Transformation to a scanner oriented coordinate frame
            cap_centres_t = cap_centres - Scanner_loc;
            cap_centres_r = (rotation_matrix * cap_centres_t')';

            %--% Elevation angles %--%
            % The vectors from the scanner to the top and bottom 'caps' of the cylinder
            cap_scan_vectors = cap_centres_r ./ sqrt(sum(cap_centres_r.^2, 2));
            
            % Their elevation angles
            cap_elev_angles = atan(cap_scan_vectors(:, num_dim) ./ sqrt(sum(cap_scan_vectors(:, 1 : num_dim - 1).^2, 2)));

            % Range of elevation angles
            elev_angle_list = min(cap_elev_angles) : laser_resolution : max(cap_elev_angles);

            %--% Azimuth angles %--%
            % The left and right outermost locations are the centres of the caps plus/minus the cylinder slice extent on the y-axis 
            cap_centres_r_rep   = repmat(cap_centres_r, [2, 1]);
            y_axis              = [0, 1, 0];
            cap_corner_points   = cap_centres_r_rep + y_axis .* [1; 1; -1; -1] * Cyl_slice_extent/2;         % TL, BL, TR, BR
            
            % The associated azimuth angles
            cap_corner_azim_angles  = atan(cap_corner_points(:, num_dim - 1) ./ cap_corner_points(:, 1));

            min_azim_angle          = min(cap_corner_azim_angles);
            max_azim_angle          = max(cap_corner_azim_angles);
            mean_azim               = (max_azim_angle + min_azim_angle) / 2;

            azim_half_amplitude     = (max_azim_angle - min_azim_angle) / 2;
            half_azim_angle_list    = (0 : laser_resolution : azim_half_amplitude) + laser_resolution;          % Adding the laser resolution acts as a round-up operator
            
            azim_angle_list         = mean_azim + [-fliplr(half_azim_angle_list), 0, half_azim_angle_list];     % The half-lists are mirrored to ensure symmetry
                        
            % The total matrices of elevation and azimuth angles
            [azim_angle_matrix, elev_angle_matrix]  = meshgrid(azim_angle_list, elev_angle_list);
            [azim_angle_list, elev_angle_list]      = deal(azim_angle_matrix(:), elev_angle_matrix(:));
                        
            % The associated vectors
            vector_matrix_r = [cos(azim_angle_list), sin(azim_angle_list), sin(elev_angle_list)];
            vector_matrix   = (rotation_matrix' * vector_matrix_r')';
            
            beam_vector_cell{s} = vector_matrix;
        end

    %% Noise-free cylinder point cloud %%                     
        point_cloud_cell_noisefree = cell(1, number_scanners);

        for s = 1 : number_scanners      
            % Directions of the laser beams
            beam_vector_matrix = beam_vector_cell{s};
                        
            % Start points of the laser beams
            Scanner_loc = Scanner_loc_cell{s};
            beam_origin_matrix = Scanner_loc .* beam_vector_matrix.^0;

            % The intersection points
            [intersection_point_matrix, cylinder_intersects, ~, incidence_angle_list] = Vector_Cylinder_Intersection(beam_origin_matrix, beam_vector_matrix, Cylinder_radius, Cylinder_centre, Cylinder_length, Cylinder_direction);

            % Only cylinder intersects with an incidence angle below the limit are of interest
            incidence_bool  = incidence_angle_list < max_incidence_angle;
            point_bool      = cylinder_intersects & incidence_bool;

            point_cloud_matrix              = intersection_point_matrix(point_bool, :);
            point_cloud_cell_noisefree{s}   = point_cloud_matrix;
        end
                
        % A structure is created containing this noise-free point cloud
        number_points_list          = cellfun(@length, point_cloud_cell_noisefree);
        Point_Cloud_Coord_Noisefree = struct('number_points_list', number_points_list, 'point_cloud_cell', {point_cloud_cell_noisefree});

        % Its uncertainty
        [sigma_radial_cell, sigma_prop_cell, ~, ~, ~] = Cylindrical_Object_Uncertainty(Cylinder_centre, Cylinder_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord_Noisefree);

        % The radial uncertainty for shifting the point in the radial direction does not come from the Gaussian power distribution, but the angular accuracy
        Angular_Accuracy_fun    = @(sigma_radial_list) angular_accuracy / beam_divergence * sigma_radial_list;
        sigma_radial_cell       = cellfun(Angular_Accuracy_fun, sigma_radial_cell, 'UniformOutput', false);

        % The distributions
        alpha                               = 1 - Confidence_interval/100;                  % The confidence interval is changed to alpha
        Point_Cloud_Distributions_Noisefree = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord_Noisefree, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Diagnostics_Distr_Generation);

    %% Simulated uncertain point cloud %%
        % They are all shifted by their noise and the distributions are resized according to any new incidence angles
        if Initial_Uncertainty == true
            % Shifted point cloud
            Point_Cloud_Coord = Gaussian_Point_Sampling(Point_Cloud_Distributions_Noisefree, Sampling_Diagnostics);

            % Its uncertainty and distributions
            [sigma_radial_cell, sigma_prop_cell, ~, ~, ~]   = Cylindrical_Object_Uncertainty(Cylinder_centre, Cylinder_direction, Scanner_loc_cell, Scanner_Parameters, Point_Cloud_Coord);
            Point_Cloud_Distributions                       = Point_Cloud_Multivariate_Normal_Generation(alpha, Point_Cloud_Coord, sigma_radial_cell, sigma_prop_cell, Scanner_loc_cell, range_bias, Diagnostics_Distr_Generation);

        % If initial uncertainty is not desired, the noisefree structures are given as output
        else
            Point_Cloud_Coord           = Point_Cloud_Coord_Noisefree;
            Point_Cloud_Distributions   = Point_Cloud_Distributions_Noisefree;
        end

    %% Uncertainty and coverage values %%
        % The distributions are projected onto the cylinder axis and cross-section
        origin                          = zeros(1, num_dim);
        [~, cylinder_vector_basis, ~]   = Vector_Based_Rotation(origin, Cylinder_direction, origin);
        
        Projected_Distributions         = Multivariate_Normal_Plane_Projection(cylinder_vector_basis, Cylinder_centre, Point_Cloud_Distributions, Diagnostics_Projection);
        
        %--% Uncertainty %--%
        % The average and maximum uncertainty of the cross-sectional ellipses and axis uncertainty are taken
        
        % Axis
        vector_sigma_cell           = Projected_Distributions.Vector.Projection.sigma;
        vector_sigma_list           = vertcat(vector_sigma_cell{:});
        
        mean_axis_uncert            = mean(vector_sigma_list);
        rel_axis_uncert             = mean_axis_uncert / Cylinder_length * 100;             % Normalised by the cylinder length
        max_axis_uncert             = max(vector_sigma_list);
        rel_max_axis_uncert         = max_axis_uncert / Cylinder_length * 100;        
        
        % Cross-section
        plane_sigmae_cell           = Projected_Distributions.Plane.Projection.sigmae;
        plane_sigmae_matrix         = vertcat(plane_sigmae_cell{:});
        proj_ellipse_sigmae_list    = sqrt(prod(plane_sigmae_matrix, 2));                   % The mean radius of an ellipse is sqrt(r_a * r_b)
        
        mean_cross_section_uncert   = mean(proj_ellipse_sigmae_list);
        rel_cross_section_uncert    = mean_cross_section_uncert / Cylinder_radius * 100;    % Normalised by the cylinder radius
        max_cross_section_uncert    = max(proj_ellipse_sigmae_list);
        rel_max_cross_sect_uncert   = max_cross_section_uncert / Cylinder_radius * 100;
        
        %--% Coverage %--%
        % Coverage is approximated numerically by determining how many of the following angles fall within at least one of the scanner's bounds
        coverage_resolution     = deg2rad(coverage_resolution);
        theta_list              = 0 : coverage_resolution : 2*pi;
        number_angles           = length(theta_list);
        
        % Coverage boolean for each scanner
        scanner_coverage_boolean = zeros(number_scanners, number_angles);
        
        for s = 1 : number_scanners
            % Scanner projected onto the cross-section
            Scanner_loc                 = Scanner_loc_cell{s};
            [~, proj_scanner_loc, ~]    = Point_to_Plane_Projection(Scanner_loc, cylinder_vector_basis, Cylinder_centre, Diagnostics_Projection);
    
            % Discretised angles w.r.t. the scanner
            theta_s     = atan2(proj_scanner_loc(2), proj_scanner_loc(1));
            phi_list    = mod(theta_list - theta_s, 2*pi);
            
            % Whether or not they fall within the scanner's coverage            
            coverage_boolean        = (phi_list < max_incidence_angle | phi_list > 2*pi - max_incidence_angle);
            scanner_coverage_boolean(s, :)  = scanner_coverage_boolean(s, :) + coverage_boolean;
        end
        
        % Resulting coverage angle
        total_coverage_boolean  = min(1, sum(scanner_coverage_boolean, 1));                 % It does not matter if one angle is covered by multiple scanners
        rel_coverage_angle      = sum(total_coverage_boolean) / number_angles * 100;        % [%]
        coverage_angle          = rel_coverage_angle / 100 * 2*pi;                          % [rad]
                
        % The length that is covered by the points
        proj_distr_height_cell  = Projected_Distributions.Vector.Projection.mu;
        proj_distr_height_list  = vertcat(proj_distr_height_cell{:});

        coverage_length     = max(proj_distr_height_list) - min(proj_distr_height_list);
        rel_coverage_length = coverage_length / Cylinder_length * 100;                              % Normalised by the cylinder length
        
        % The results are saved in a structure and file
        number_beams_total      = sum(number_points_list);
        Uncertainty_Coverage    = struct('number_beams_total', number_beams_total, 'mean_cross_section_uncert', mean_cross_section_uncert, 'rel_cross_section_uncert', rel_cross_section_uncert, 'max_cross_section_uncert', max_cross_section_uncert, 'rel_max_cross_section_uncert', rel_max_cross_sect_uncert, ... 
                                         'mean_axis_uncert', mean_axis_uncert, 'rel_axis_uncert', rel_axis_uncert, 'max_axis_uncert', max_axis_uncert, 'rel_max_axis_uncert', rel_max_axis_uncert, 'coverage_angle', coverage_angle, 'rel_coverage_angle', rel_coverage_angle, 'coverage_length', coverage_length, 'rel_coverage_length', rel_coverage_length);
        
        save('Point_Cloud_Uncertainty_and_Coverage.mat', 'Uncertainty_Coverage');
                                     
        % Printed messages
        if Print == true
            disp('-----------------------------------------------------------------------')
            fprintf('%g synthetic points were generated \n',            number_beams_total);
            disp('The average uncertainty w.r.t. the:');
            fprintf('   Cross-section: %.3g mm, %.3g%% of radius \n',  mean_cross_section_uncert * 1e3, rel_cross_section_uncert);         % Note conversion to mm
            fprintf('   Axis:          %.3g mm, %.3g%% of length \n',  mean_axis_uncert * 1e3, rel_axis_uncert);                           % Note conversion to mm
            disp('The maximum uncertainty w.r.t. the:');
            fprintf('   Cross-section: %.3g mm, %.3g%% of radius \n',  max_cross_section_uncert * 1e3, rel_max_cross_sect_uncert);         % Note conversion to mm
            fprintf('   Axis:          %.3g mm, %.3g%% of length \n',  max_axis_uncert * 1e3, rel_max_axis_uncert);                        % Note conversion to mm
            disp('The coverage w.r.t. the:');
            fprintf('   Cross-section: %.3g deg, %.3g%% of full circle \n',           rad2deg(coverage_angle), rel_coverage_angle);        % Note conversion to deg
            fprintf('   Axis:          %.3g m, %.3g%% of length \n',   coverage_length, rel_coverage_length);
            disp('-----------------------------------------------------------------------')           
        end
        
    %% Plot %% 
        if Diagnostics == true
            %--% Plot %--%
            % Number of standard deviations shown by the distributions
            number_STD      = 2;

            % The colours used for the beams and uncertainties of each scanner
            scanner_cmap    = cbrewer('qual', 'Set1', max(3, number_scanners));
            
            % The distributions
            distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;
            distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
            distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;

            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])    
            
            hold on
            grid off
                
            % The cylinder surface
            number_coord = 1e3;
            [cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, ~] = Cylinder_Surface_Generator(Cylinder_radius, Cylinder_length, Cylinder_centre, Cylinder_direction, number_coord);
            surf(cylinder_coord_x, cylinder_coord_y, cylinder_coord_z, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.25, 'LineWidth', 2, 'DisplayName', 'Cylinder');

            % Its caps
            pl_top = plot3(cylinder_coord_x(1, :), cylinder_coord_y(1, :), cylinder_coord_z(1, :), 'LineWidth', 2, 'color', 'k');
            pl_bot = plot3(cylinder_coord_x(2, :), cylinder_coord_y(2, :), cylinder_coord_z(2, :), 'LineWidth', 2, 'color', 'k');            
            pl_top.HandleVisibility = 'Off';
            pl_bot.HandleVisibility = 'Off';
            
            % Each scanner's beams, points and distributions
            point_cloud_cell        = Point_Cloud_Coord.point_cloud_cell;

            number_cumulative_beams = [0, cumsum(number_points_list)];

            for s = 1 : number_scanners
                % This scanner's data
                Scanner_loc     = Scanner_loc_cell{s};
                scanner_colour  = scanner_cmap(s, :);

                ind_first       = number_cumulative_beams(s) + 1;
                ind_last        = number_cumulative_beams(s + 1);
                
                point_cloud_matrix  = point_cloud_cell{s};

                % Locations of the points
                scanner_range = norm(Scanner_loc - Cylinder_centre);
                scatter3(point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'MarkerFaceColor', scanner_colour, 'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Point cloud, scanner %i (R = %i m)', s, scanner_range));
                
                for ind = ind_first : ind_last
                    % Distribution's properties
                    distribution_mu     = distribution_mu_cell{ind};
                    distribution_radii  = number_STD * distribution_sigmae_cell{ind};
                    distribution_axes   = distribution_axes_cell{ind};

                    % Laser beam
                    [mu_x, mu_y, mu_z]  = Column_Deal(distribution_mu);

                    laser_beam_string   = sprintf('Laser beams, scanner %i', s);
                    pl_laser            = plot3([Scanner_loc(1), mu_x], [Scanner_loc(2), mu_y], [Scanner_loc(3), mu_z], 'color', scanner_colour, 'DisplayName', laser_beam_string);

                    % Distribution
                    [ellipsoid_coord_matrix, number_coord] = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);

                    x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);

                    surf_el = surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', scanner_colour, 'FaceAlpha', 0.25, 'DisplayName', sprintf('Fuzzy cloud (%i %s), scanner %i', number_STD, '\sigma', s));

                    if ind > ind_first
                        pl_laser.HandleVisibility   = 'Off';
                        surf_el.HandleVisibility    = 'Off';
                    end
                end
            end
                        
            % Axes
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('z [m]')

            axis equal

            point_cloud_matrix_tot = vertcat(point_cloud_cell{:});

            x_lim = Cylinder_radius * [-1, 1] + [min(point_cloud_matrix_tot(:, 1)), max(point_cloud_matrix_tot(:, 1))];
            xlim(x_lim);
            y_lim = Cylinder_radius * [-1, 1] + [min(point_cloud_matrix_tot(:, 2)), max(point_cloud_matrix_tot(:, 2))];
            ylim(y_lim);            
            z_lim = Cylinder_radius * [-1, 1] + [min(point_cloud_matrix_tot(:, 3)), max(point_cloud_matrix_tot(:, 3))];
            zlim(z_lim);     
            
            % Viewing angle
            view(120, 25)

            % Legend
            legend('show', 'location', 'eastoutside');

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off

            disp('The script will continue and the plot will close when a key is pressed')
            pause()

            close(1)            
        end
        
end