% The top and bottom of the cylinder (and thus the length) are determined by integrating the point cloud probability distributions onto the cylinder axis.
% These PDFs are then weighted by their repulsion vectors, and the maximum likelihoods correspond to the most likely edge locations

function [Cylinder_length, Cylinder_top_loc, Cylinder_bot_loc] = Cylinder_Length_Estimation(Cylinder_centre, Cylinder_direction, Fitting_Parameters, Point_Cloud_Distributions, Output_Decisions)

    %% Inputs %%
        % Point cloud distributions
        distribution_axes_cell      = Point_Cloud_Distributions.distribution_axes_cell;
        distribution_sigmae_cell    = Point_Cloud_Distributions.distribution_sigmae_cell;
        distribution_Sigma_cell     = Point_Cloud_Distributions.distribution_Sigma_cell;
        distribution_mu_cell        = Point_Cloud_Distributions.distribution_mu_cell;

        % Fitting parameters
        point_weights_list          = Fitting_Parameters.point_weights_list;
        alignment_threshold         = Fitting_Parameters.alignment_threshold;

        % Output decisions
        Plot                        = Output_Decisions.Plot;
        Print                       = Output_Decisions.Print;

    %% Manual inputs %%
        % Intermediate edge detection and plane fitting plots can be created for diagnostic purposes
        Diagnostics     = false;            % [true, false]
        Edge_Plot       = false;            % [true, false]
        Fitting_Plot    = false;            % [true, false]

    %% Edge detection %%
        % Top and bottom edge points are detected using a 3D triangular mesh of the point cloud
        [edge_classes, number_edge_classes, point_cloud_edge_bool_f, ~, number_edge_points_list, surface_triangle_mesh] = Cylinder_Mesh_Edge_Points(alignment_threshold, Cylinder_direction, Cylinder_centre, Point_Cloud_Distributions, Diagnostics, Edge_Plot);
        
        % The point cloud distributions structure is now defined for each edge
        edge_PC_distributions_cell = cell(1, number_edge_classes);
    
        for e = 1 : number_edge_classes
            % This edge's structure
            edge_bool                       = point_cloud_edge_bool_f(e, :);
            Edge_Point_Cloud_Distributions  = struct('distribution_mu_cell', {distribution_mu_cell(edge_bool)}, 'distribution_sigmae_cell', {distribution_sigmae_cell(edge_bool)}, 'distribution_Sigma_cell', {distribution_Sigma_cell(edge_bool)}, 'distribution_axes_cell', {distribution_axes_cell(edge_bool)}, 'number_distributions', number_edge_points_list(e));
            edge_PC_distributions_cell{e}   = Edge_Point_Cloud_Distributions;
        end

    %% Heights of the top and bottom planes %%
        % The plane is defined as orthogonal to the given cylinder direction, so only the height needs to be determined
        Plane_Height = struct();
    
        for e = 1 : number_edge_classes
            % This edge's data
            number_edge_points  = number_edge_points_list(e);
            edge_class          = edge_classes{e};

            % If there are no points in this edge class, the lowest/highest mu is used
            if number_edge_points == 0
                point_cloud_matrix  = vertcat(distribution_mu_cell{:});
                [~, delta_list, ~]  = Point_to_Vector_Projection(point_cloud_matrix, Cylinder_direction, Cylinder_centre);

                if strcmp(edge_class, 'Top')
                    Plane_Height.(edge_class) = max(delta_list);
                elseif strcmp(edge_class, 'Bottom')
                    Plane_Height.(edge_class) = min(delta_list);
                end
            else
                edge_bool                       = point_cloud_edge_bool_f(e, :);
                Edge_Point_Cloud_Distributions  = edge_PC_distributions_cell{e};
    
                % The point weights of these distributions
                edge_point_weights_list = point_weights_list(edge_bool);
                edge_point_weights_list = edge_point_weights_list / mean(edge_point_weights_list);
    
                Edge_Fitting_Parameters                     = Fitting_Parameters;
                Edge_Fitting_Parameters.point_weights_list  = edge_point_weights_list;
    
                % The optimal plane height with respect to the origin
                Fixed_Normal_Vector             = true;             
                [~, ~, ~, plane_height, ~, ~]   = Fuzzy_Infinite_Plane_Fitting(Edge_Point_Cloud_Distributions, Edge_Fitting_Parameters, Fixed_Normal_Vector, Cylinder_direction, Diagnostics, Fitting_Plot);
                Plane_Height.(edge_class)       = plane_height;
            end
        end

    %% Cylinder length %%
        % Cylinder top and bottom heights and locations
        Cylinder_Edges = struct();

        for e = 1 : number_edge_classes
            % The edge class (top or bottom)
            edge_class      = edge_classes{e};

            % Plane location
            plane_height    = Plane_Height.(edge_class);
            plane_location  = plane_height * Cylinder_direction;

            % The cylinder edge locations is its projection onto the cylinder axis
            [cylinder_edge_loc, cylinder_edge_height, ~]    = Point_to_Vector_Projection(plane_location, Cylinder_direction, Cylinder_centre);
            Cylinder_Edges.(edge_class).Location            = cylinder_edge_loc;
            Cylinder_Edges.(edge_class).Height              = cylinder_edge_height;
        end
        
        % The top and bottom locations
        Cylinder_top_loc        = Cylinder_Edges.Top.Location;
        Cylinder_bot_loc        = Cylinder_Edges.Bottom.Location;

        % The length
        Cylinder_top_height     = Cylinder_Edges.Top.Height;
        Cylinder_bottom_height  = Cylinder_Edges.Bottom.Height;
        Cylinder_length         = Cylinder_top_height - Cylinder_bottom_height;
        
    %% Fitting assessment %%
        % Printed messages
        if Print == true
            fprintf('The cylinder has been fitted with length %.3g m \n', Cylinder_length);
            fprintf('And top:    [%.3g, %.3g, %.3g] \n', Cylinder_top_loc(1), Cylinder_top_loc(2), Cylinder_top_loc(3));
            fprintf('And bottom: [%.3g, %.3g, %.3g] \n', Cylinder_bot_loc(1), Cylinder_bot_loc(2), Cylinder_bot_loc(3));
        end

        % Plot 
        if Plot == true
            % Edge colour map
            edge_cmap = cbrewer('qual', 'Set1', max(number_edge_classes, 3));
            edge_cmap = max(edge_cmap, 0);
            edge_cmap = min(edge_cmap, 1);

            % Point cloud matrix
            point_cloud_matrix = vertcat(distribution_mu_cell{:});

            % Plot settings
            number_coord    = 1e2;
            m_STD           = 1;
    
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])
    
            for e = 1 : number_edge_classes
                % This edge's data
                edge_class          = edge_classes{e};
                edge_indices        = find(point_cloud_edge_bool_f(e, :));
                edge_colour         = edge_cmap(e, :);
                number_edge_points  = number_edge_points_list(e);

                %--% Subplot %--%
                subplot(1, number_edge_classes, e)
                hold on
                grid on            
                
                % The surface triangles
                trimesh(surface_triangle_mesh, point_cloud_matrix(:, 1), point_cloud_matrix(:, 2), point_cloud_matrix(:, 3), 'LineWidth', 1, 'EdgeColor', 'k', 'FaceColor', [1, 1, 1], 'FaceAlpha', 0.5, 'DisplayName', 'Surface mesh');
    
                % The points belonging to the surface facets       
                surface_ind            = surface_triangle_mesh(:);
                surface_coordinates    = point_cloud_matrix(surface_ind, :);
                scatter3(surface_coordinates(:, 1), surface_coordinates(:, 2), surface_coordinates(:, 3), 'filled', 'g', 'DisplayName', 'Surface points');
        
                % Edge distributions
                scatter3(point_cloud_matrix(edge_indices, 1), point_cloud_matrix(edge_indices, 2), point_cloud_matrix(edge_indices, 3), 50, 'filled', 'MarkerFaceColor', edge_colour, 'DisplayName', sprintf('%s_{%s}', '\mu', edge_class));
            
                for i = 1 : number_edge_points
                    % Distribution properties
                    ind                 = edge_indices(i);
                    distribution_mu     = distribution_mu_cell{ind};
                    distribution_sigmae = distribution_sigmae_cell{ind};
                    distribution_radii  = distribution_sigmae * m_STD;
                    distribution_axes   = distribution_axes_cell{ind};
    
                    % Its coordinates at m sigma
                    [distr_coord_matrix, number_coord]  = Ellipsoid_Coordinate_Generator(distribution_mu, distribution_radii, distribution_axes, number_coord);
    
                    x_distribution  = reshape(distr_coord_matrix(:, 1), sqrt(number_coord) * [1, 1]);
                    y_distribution  = reshape(distr_coord_matrix(:, 2), sqrt(number_coord) * [1, 1]);
                    z_distribution  = reshape(distr_coord_matrix(:, 3), sqrt(number_coord) * [1, 1]);
    
                    % Its surface patch
                    surf_distr_name = sprintf('Distribution, %.2g %s', m_STD, '\sigma');
                    distr_surf = surf(x_distribution, y_distribution, z_distribution, 'EdgeColor', 'none', 'FaceColor', edge_colour, 'FaceAlpha', 0.25, 'DisplayName', surf_distr_name);
    
                    if i > 1
                        distr_surf.HandleVisibility = 'Off';
                    end
                end

                % The fitted plane
                cylinder_edge_loc   = Cylinder_Edges.(edge_class).Location;
                plane_corner_matrix = Plane_Corner_Points(Cylinder_direction, cylinder_edge_loc, point_cloud_matrix);
                patch(plane_corner_matrix(:, 1), plane_corner_matrix(:, 2), plane_corner_matrix(:, 3), edge_colour, 'FaceAlpha', 0.25, 'DisplayName', sprintf('%s plane', edge_class));   

                % Axes
                xlabel('x [m]');
                ylabel('y [m]');
                zlabel('z [m]');
    
                axis equal
    
                % Viewing angle
                view(45, 45);
    
                % Legend
                legend('show', 'location', 'northoutside');
    
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off    
            end
    
            % Pause message
            disp('The cylinder length has been determined. The figures will close and script end upon a key-press');
            pause();
    
            close(1);
        end
end