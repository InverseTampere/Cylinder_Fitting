% This script computes the radial and range uncertainty of laser beams hitting a circular object, where the circle has centre point [circle_x, circle_y] and radius circle_r
% The beams originate from [x_scanner, y_scanner] to hits [x_list, y_list] with the given beam divergence
% Uncertainty in the range estimate due to the device can also be added

% The incidence angle (and thus range uncertainty) can grow to infinity for grazing hits if the circle geometry is estimated
% As such, for points for which the maximum incidence angle is exceeded the uncertainty is substituted
% Max_incidence_angle can bet set to >pi/2 if this is undesired

function [sigma_radial_cell, sigma_range_cell, incidence_angle_cell, beam_range_cell] = Circular_Object_Uncertainty(circle_x, circle_y, circle_r, Scanner_loc_cell, x_points_cell, y_points_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication)

    %% The incidence angle of each beam %%
        number_scanners = length(Scanner_loc_cell);

        incidence_angle_cell = cell(1, number_scanners);
        
        for s = 1 : number_scanners
            % The data for this scanner 
            x_points_list   = x_points_cell{s};
            y_points_list   = y_points_cell{s};
            
            % If the point cloud was duplicated, only the first half is used and the uncertainty itself is duplicated
            if Point_Duplication == true
                number_points   = length(x_points_list);
                x_points_list   = x_points_list(1 : number_points / 2);
                y_points_list   = y_points_list(1 : number_points / 2);
            end
            
            x_scanner       = Scanner_loc_cell{s}(1);
            y_scanner       = Scanner_loc_cell{s}(2);
            
            % The beam vectors
            vector_x_list       = x_points_list - x_scanner;
            vector_y_list       = y_points_list - y_scanner;
            
            % The intersection points of the laser beams with the curved surface
            number_beams        = length(x_points_list);
            x_start_list        = x_scanner * ones(1, number_beams);
            y_start_list        = y_scanner * ones(1, number_beams);

            [x_i_hits, y_i_hits, ~, ~] = Vector_Circle_Intersection(x_start_list, y_start_list, vector_x_list, vector_y_list, circle_x, circle_y, circle_r);

            % The incidence angles between the beams and the curved surface
            incidence_angle_list    = Curved_Object_Incidence_Angle(x_i_hits, y_i_hits, vector_x_list, vector_y_list, circle_x, circle_y);

            % To ensure that excessively large uncertainties are not present, the incidence angle is capped
            incidence_angle_list    = min(incidence_angle_list, max_incidence_angle);
            incidence_angle_cell{s} = incidence_angle_list;
        end
        
    %% The uncertainty of the beams for each scanner %%        
        [sigma_radial_cell, sigma_range_cell, beam_range_cell] = Laser_Beam_Uncertainty(x_points_cell, y_points_cell, incidence_angle_cell, Scanner_loc_cell, beam_divergence, sigma_range_device, max_incidence_angle, Point_Duplication);

    %% The results are duplicated if need be %%
        if Point_Duplication == true
            duplicate_fun           = @(x) [x; x];
            
            sigma_radial_cell       = cellfun(duplicate_fun, sigma_radial_cell, 'UniformOutput', false);
            sigma_range_cell        = cellfun(duplicate_fun, sigma_range_cell, 'UniformOutput', false);
            incidence_angle_cell    = cellfun(duplicate_fun, incidence_angle_cell, 'UniformOutput', false);
            beam_range_cell         = cellfun(duplicate_fun, beam_range_cell, 'UniformOutput', false);
        end
end