% This script computes the objective value of the cross-section fit using the Mahalanobis distance between the point cloud and the 3D circle

function Objective_value = Cylinder_Cross_Section_Objective(inf_cylinder_geometry_n, Geometry_Indices, geometry_LB, geometry_UB, Objective_value_init, Fitting_Parameters, Point_Cloud_Distributions_c, Scanning_Parameters_c, Statistical_Values)
        
    %% Structure inputs %%
        % Geometry indices
        centre_ind              = Geometry_Indices.centre;
        radius_ind              = Geometry_Indices.radius;
        azim_ind                = Geometry_Indices.azim_angle;
        elev_ind                = Geometry_Indices.elev_angle;

        % Fitting parameters
        distance_moment         = Fitting_Parameters.distance_moment;
        Point_weights_list      = Fitting_Parameters.point_weights_list;
        Distance_Computation    = Fitting_Parameters.Distance_Computation;

        % Scanner locations
        Scanner_loc_cell_c      = Scanning_Parameters_c.Scanner_loc_cell;

    %% Manual inputs %%
        % It is recommend that outputs are set to false
        Diagnostics             = false;        % [true, false]
   
    %% Geometry parameters (optimisation variables) %%
        % Un-normalising the geometry vector
        inf_cylinder_geometry = inf_cylinder_geometry_n .* (geometry_UB - geometry_LB) + geometry_LB;

        % Circular cross section
        circle_centre   = inf_cylinder_geometry(centre_ind);               % Expressed in the cylinder-aligned coordinate frame (i.e. projection on the plane orthogonal to the cylinder axis crossing the origin)
        circle_radius   = inf_cylinder_geometry(radius_ind);    
        
        % Cylinder direction from the spherical angles
        azim_angle              = inf_cylinder_geometry(azim_ind);
        azim_angle              = azim_angle - sign(azim_angle)*2*pi;

        elev_angle              = inf_cylinder_geometry(elev_ind);
        
        if elev_angle < 0
            elev_angle          = abs(elev_angle);
        elseif elev_angle > pi/2
            elev_angle          = pi - elev_angle;
        end

        [cylinder_dir, ~, ~]    = Vector_Spherical_Angle_Conversion([], azim_angle, elev_angle);

        % The circle centre is converted to the cylinder centre
        num_dim = length(cylinder_dir);
        origin  = zeros(1, num_dim);                                                            % The origin is used as the centroid as the point cloud is already centered
        
        circle_height                   = 0;                                                    % Note that the height above the cross-section is irrelevant when the problem is projected onto it
        [~, ~, ~, cylinder_centre_c]    = Circle_Cylinder_Centre_Conversion(circle_centre, circle_height, [], cylinder_dir, origin, Scanner_loc_cell_c);

    %% The normalised expected Mahalanobis distance of the point cloud given this cross-section %%
        % The weighted expected (squared) Mahalanobis distance
        if strcmp(Distance_Computation, 'Line_Approx')
            [~, expected_Mahal_distance_list] = Expected_Mahal_Distance_Cylinder_Line_Approx(cylinder_centre_c, circle_radius, cylinder_dir, distance_moment, Point_Cloud_Distributions_c, Scanner_loc_cell_c, Diagnostics);
        elseif strcmp(Distance_Computation, 'Numerical')
            [~, expected_Mahal_distance_list] = Expected_Mahalanobis_Distance_Cylinder_Numerical(cylinder_centre_c, circle_radius, cylinder_dir, distance_moment, Point_Cloud_Distributions_c, Statistical_Values, Diagnostics);
        end
                
        weighted_expected_Mahal_distance_list = Point_weights_list .* expected_Mahal_distance_list;

        % Normalised by the value for the initial geometry to form the objective function
        Objective_value_list = weighted_expected_Mahal_distance_list / Objective_value_init;
        Objective_value      = mean(Objective_value_list);

    %% Printed messages %%
        % For diagnostic purposes the objective value, gradient and geometry are printed
        if Diagnostics == true
            disp('---');
            fprintf('O = %.3g, E[M^2] = %.3g \n', Objective_value, mean(weighted_expected_Mahal_distance_list));
            fprintf('c = [%.3g, %.3g, %.3g] m \n', cylinder_centre_c);
            fprintf('v = [%.3g, %.3g, %.3g] m \n', cylinder_dir);
            fprintf('r = %.3g m \n', circle_radius);       
        end
end
        