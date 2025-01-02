% This script optimises the circle fit through a combination of the normalised likelihood and the expected distance of the given moment between the point cloud and the circle
% The normalised likelihood is used, as it is independent of sigma and laser beamwidth and ranges between 0 and 1
% As optimisation is performed through minimisation, the complement of the normalised likelihood is used (i.e. Lc = 1 - L)

% The likelihood and distance objectives are balanced according to objective_balance, which is the factor by which the distance objective is multiplied

% If weighting is true, weights are the inverse of the expected distance of each point

function Objective_value = Circle_Fitting_Optimiser(objective_balance, distance_moment, Tree_stem_geometry_n, geometry_lb, geometry_ub, x_points_cell, y_points_cell, Expected_distance_cell_init, sigma_radial_init_cell, sigma_range_init_cell, laser_radius_init_cell, Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Weighting, central_theta_list, delta_inner_list, delta_outer_list)

    %% Tree stem parameters (optimisation variables) %%
        Tree_stem_geometry  = Tree_stem_geometry_n .* (geometry_ub - geometry_lb) + geometry_lb;

        Tree_centre_x       = Tree_stem_geometry(1);    % Centre of the tree stem
        Tree_centre_y       = Tree_stem_geometry(2);
        Tree_radius         = Tree_stem_geometry(3);    % The tree stem's radius
        
    %% The average expected distance of this point cloud to the circle %%
        if objective_balance ~= 0                       % If this objective is not accounted for, computing it is not useful
            Plot    = false;
            Print   = false;
            
            [Average_expected_distance, ~]  = Expected_Circle_Distance(distance_moment, Scanner_loc_cell, x_points_cell, y_points_cell, sigma_radial_init_cell, sigma_range_init_cell, range_bias, laser_radius_init_cell, Tree_centre_x, Tree_centre_y, Tree_radius, Point_Mirroring, Point_Duplication, Plot, Print);
            
            % Normalised to better fit the optimiser       
            Expected_distance_list_init     = horzcat(Expected_distance_cell_init{:});
            Expected_distance_init          = mean(Expected_distance_list_init);

            Norm_expected_distance          = Average_expected_distance / Expected_distance_init;
        else
            Norm_expected_distance          = 0;
        end
                
    %% The normalised likelihood of this circle fitting the given point cloud %%
        Diagnostics = false;        % Displaying things during optimisation quickly becomes overwhelming
        [~, ~, Avg_L_N_max, L_N_max_cell, ~, ~] = Circle_Fit_Likelihood(x_points_cell, y_points_cell, sigma_radial_init_cell, sigma_range_init_cell, laser_radius_init_cell, Tree_centre_x, Tree_centre_y, Tree_radius, Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Diagnostics);
        
        %--% The normalised likelihood is weighted (if desired) %--%
        if Weighting == true                        
            % The weights are the inverse of the squared expected distance
            Weights_list                = 1./Expected_distance_list_init;
            
            % The weights are normalised to have a max of 1
            Weights_list                = Weights_list / max(Weights_list);
            
            % The normalised likelihoods are multiplied by their respective weights
            Norm_likelihood_list        = horzcat(L_N_max_cell{:});
            Weighted_Likelihood_N_list  = Norm_likelihood_list .* Weights_list;
            
            Avg_L_N_max         = mean(Weighted_Likelihood_N_list);
        end

    %% The objective function %%
        [Geometry_Function_Handle, Geometry_variable_bounds] = Shape_Basis_Function();
    
        Objective_value = Objective_Function(Avg_L_N_max, Norm_expected_distance, objective_balance, Geometry_Function_Handle, Tree_stem_geometry, Geometry_variable_bounds, central_theta_list, delta_inner_list, delta_outer_list, Scanner_loc_cell);

end
