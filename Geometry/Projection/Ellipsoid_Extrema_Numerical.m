% This script computes the extrema of a 3D ellipsoid, for the desired dimensions (i.e. its minimum and maximum x, y and z coordinates)
% Note that the extrema for dimensions outside of Dimensions are not computed and are given as NaN

% It does this through optimisation, as it can be seen as a bounded 2-parameter optimisation problem (spherical parametrisation of the ellipsoid)

function [extrema_matrix, extrema_loc_matrix] = Ellipsoid_Extrema_Numerical(Dimensions, ellipsoid_centre, ellipsoid_radii, ellipsoid_axes)

    %% Manual inputs %%
        Diagnostics     = false;             % [true, false] Shows the ellipsoid and its extrema

    %% Finding the extrema %%
        % The optimisation options
        Options         = optimoptions('fmincon', 'optimalitytolerance', 1e-4, 'steptolerance', 1e-4, 'Display', 'Off');
            
        % The optimisation inputs
        Optim_LB    = [0, 0];
        Optim_UB    = [pi, 2*pi];
        Init_Angles = (Optim_LB + Optim_UB)/2;     
        
        % The objective functions handle
        Objective_fun_handle    = @(Angles) Ellipsoid_Delta_fun(Angles);
        
        % The extrema in each desired dimension  
        num_dim = length(ellipsoid_centre);

        extrema_matrix      = NaN(num_dim, 2);              % [min, max]
        extrema_loc_matrix  = NaN(num_dim, num_dim, 2);     % Locations of [min, max] in 3rd dimension
        
        for dim = Dimensions
            % Optimisation
            [Extremum_Angles, ~, ~] = fmincon(Objective_fun_handle, Init_Angles, [], [], [], [], Optim_LB, Optim_UB, [], Options);
            
            % The resulting extrema
            % Note that an ellipsoid is symmetrical w.r.t. its centroid
            ellipsoid_extremum_delta_coord  = Ellipsoid_Delta_Coord_fun(Extremum_Angles, ellipsoid_radii, ellipsoid_axes);
            extrema_loc_matrix(dim, :, :)   = ellipsoid_centre + cat(3, -1, 1) .* ellipsoid_extremum_delta_coord;
            
            extrema_matrix(dim, :)          = extrema_loc_matrix(dim, dim, :);
        end
    
    %% Extrema function %%
        function delta_dimension_value = Ellipsoid_Delta_fun(Angles)
            % Resulting location on the ellipsoid
            ellipsoid_delta_coord = Ellipsoid_Delta_Coord_fun(Angles, ellipsoid_radii, ellipsoid_axes);
            
            % The desired dimension value (w.r.t. ellipsoid centre)
            delta_dimension_value = ellipsoid_delta_coord(dim);
        end
    
    %% Ellipsoid coordinate function %%
        function ellipsoid_delta_coord = Ellipsoid_Delta_Coord_fun(Angles, ellipsoid_radii, ellipsoid_axes)
            % The angles
            alpha   = Angles(1);
            beta    = Angles(2);
            
            % The resulting location on the ellipsoid (w.r.t. ellipsoid centre)
            unit_sphere_coord       = [sin(alpha)*cos(beta); sin(alpha)*sin(beta); cos(alpha)];
            ellipsoid_unrot_coord   = ellipsoid_radii' .* unit_sphere_coord;
            
            ellipsoid_delta_coord   = (ellipsoid_axes \ ellipsoid_unrot_coord)';
        end
    
    %% Plot %%
        if Diagnostics == true
            % Ellipsoid coordinates
            m_ellipsoid = 1e3;              % The number of points in the ellipse [-]
            [ellipsoid_coord_matrix, m_ellipsoid] = Ellipsoid_Coordinate_Generator(ellipsoid_centre, ellipsoid_radii, ellipsoid_axes, m_ellipsoid);

            % Text used for plotting
            dimension_names = {'i', 'j', 'k'};          % Generic indications for the dimensions are given
            extrema_names   = {'Minima', 'Maxima'};
            num_extrema     = length(extrema_names);

            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Ellipsoid
            x_ellipsoid = reshape(ellipsoid_coord_matrix(:, 1), sqrt(m_ellipsoid) * [1, 1]);
            y_ellipsoid = reshape(ellipsoid_coord_matrix(:, 2), sqrt(m_ellipsoid) * [1, 1]);
            z_ellipsoid = reshape(ellipsoid_coord_matrix(:, 3), sqrt(m_ellipsoid) * [1, 1]);

            surf(x_ellipsoid, y_ellipsoid, z_ellipsoid, 'EdgeColor', 'none', 'FaceColor', 'g', 'FaceAlpha', 0.1, 'DisplayName', 'Ellipsoid');

            % Extrema
            for e = 1 : num_extrema
                % The name of this extremum
                extremum_name = extrema_names{e};

                for d = 1 : num_dim
                    % This dimension
                    dimension = dimension_names{d};

                    legend_string = sprintf('%s %s', extremum_name, dimension);
                    scatter3(extrema_loc_matrix(d, 1, e), extrema_loc_matrix(d, 2, e), extrema_loc_matrix(d, 3, e), 'filled', 'DisplayName', legend_string); 
                end
            end

            % Axes
            xlabel('x [m]');
            ylabel('y [m]');
            zlabel('z [m]');

            max_list    = max(ellipsoid_coord_matrix, [], 1);
            min_list    = min(ellipsoid_coord_matrix, [], 1);
            ampl_list   = max_list - min_list;

            xlim(0.2*ampl_list(1)*[-1, 1] + [min_list(1), max_list(1)]);
            ylim(0.2*ampl_list(2)*[-1, 1] + [min_list(2), max_list(2)]);
            zlim(0.2*ampl_list(3)*[-1, 1] + [min_list(3), max_list(3)]);

            axis equal      % Note that this changes the axis limits again

            view(45, 45);

            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off    
            
            % Pause
            disp('The script will continue and the plot will close upon a button-press.');
            pause();
            
            close(1);
        end
end