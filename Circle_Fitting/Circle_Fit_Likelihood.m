% This script computes the normalised likelihood of the given circle, given an uncertain point cloud

% The normalised likelihood is used, as it is independent of sigma and laser beamwidth and ranges between 0 and 1
% The precision to which this value is found for each point is at least L_N_precision
% Additionally, as the distance between the circle and point are lower for points inside the circle, they can be mirrored by setting Point_Mirroring to true

function [Avg_L_N_total, L_N_total_cell, Avg_L_N_max, L_N_max_cell, Avg_L_N_min, L_N_min_cell] = Circle_Fit_Likelihood(x_points_cell, y_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Tree_centre_x, Tree_centre_y, Tree_radius, Scanner_loc_cell, range_bias, Point_Duplication, Point_Mirroring, Diagnostics)

    %% Point cloud preprocessing %%       
        %--% Range bias compensation %--%
        [x_points_cell, y_points_cell] = Range_Bias_Compensation(x_points_cell, y_points_cell, range_bias, Scanner_loc_cell);

        %--% Point cloud mirroring %--%
        if Point_Mirroring == true           
            number_scanners     = length(x_points_cell);

            x_int_points_cell   = cell(1, number_scanners);
            y_int_points_cell   = cell(1, number_scanners);

            for s = 1 : number_scanners
                % This scanner's location and point cloud
                Scanner_loc_x   = Scanner_loc_cell{s}(1);
                Scanner_loc_y   = Scanner_loc_cell{s}(2);

                x_points_list   = x_points_cell{s};
                y_points_list   = y_points_cell{s};

                number_beams    = length(x_points_list);

                %--% Intersection points %--%
                % For beams that intersect the circle, the intersection point is used  
                x_start_list        = Scanner_loc_x * ones(number_beams, 1);
                y_start_list        = Scanner_loc_y * ones(number_beams, 1);

                vector_beams_x_list = x_points_list - x_start_list;
                vector_beams_y_list = y_points_list - y_start_list;

                [x_i_list, y_i_list, ~, intersection_list]  = Vector_Circle_Intersection(x_start_list, y_start_list, vector_beams_x_list, vector_beams_y_list, Tree_centre_x, Tree_centre_y, Tree_radius);

                % If no intersection point is used, the closest point on the circle to the point is used
                missed_beams = intersection_list == false;

                % The intersection points are the closest point on the circle to the point cloud, i.e. from the points to the tree stem centre
                vector_point_tree_x_list    = x_points_list - Tree_centre_x;
                vector_point_tree_y_list    = y_points_list - Tree_centre_y;

                [x_orth_i_list, y_orth_i_list, delta_orth_list, ~]  = Vector_Circle_Intersection(x_points_list, y_points_list, vector_point_tree_x_list, vector_point_tree_y_list, Tree_centre_x, Tree_centre_y, Tree_radius);

                % These intersection points replace the incorrect values
                x_i_list(missed_beams) = x_orth_i_list(missed_beams);
                y_i_list(missed_beams) = y_orth_i_list(missed_beams);

                % If the points were duplicated, the intersection points must be duplicated as well
                if Point_Duplication == true
                    % The first half belongs to the original points
                    x_i_first_half_list = x_i_list(1 : number_beams/2);
                    y_i_first_half_list = y_i_list(1 : number_beams/2);

                    % The intersection points are duplicated
                    [x_i_duplicated_cell, y_i_duplicated_cell] = Point_Cloud_Duplication(Tree_centre_x, Tree_centre_y, {x_i_first_half_list}, {y_i_first_half_list}, 1);            
                    x_i_list = x_i_duplicated_cell{1};
                    y_i_list = y_i_duplicated_cell{1};

                    % The distances are the same
                    delta_orth_list(number_beams/2 + 1 : end) = delta_orth_list(1 : number_beams/2);
                end

                x_int_points_cell{s} = x_i_list;
                y_int_points_cell{s} = y_i_list;

                %--% The points are mirrored  %--%
                % It is ensured that all points lie on the same side of the circle to avoid the asymmetrical distances affecting the fit by mirroring them
                % The points are moved from inside to outside
                % When moved inside the average distance is lower, and thus the area of non-zero likelihood is greater as well
                outlier_points = delta_orth_list < 0;

                x_points_list(outlier_points) = 2*x_orth_i_list(outlier_points) - x_points_list(outlier_points);
                y_points_list(outlier_points) = 2*y_orth_i_list(outlier_points) - y_points_list(outlier_points);    

                x_points_cell{s} = x_points_list;
                y_points_cell{s} = y_points_list;
            end
        end
        
    %% The likelihood %%
        [Avg_L_N_total, L_N_total_cell, Avg_L_N_max, L_N_max_cell, Avg_L_N_min, L_N_min_cell] = Total_Circle_Likelihood(x_points_cell, y_points_cell, sigma_radial_cell, sigma_range_cell, laser_radius_cell, Tree_centre_x, Tree_centre_y, Tree_radius, Scanner_loc_cell, Diagnostics, Point_Mirroring);
    
end


    