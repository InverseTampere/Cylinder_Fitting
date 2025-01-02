% The viewing angles from which the laser scanners can see the tree for the given geometry and resolution
% Additionally, the limiting angles are given

function [theta_cell, theta_LB_cell, theta_UB_cell] = Scanning_Angles(Tree_centre, Tree_radius, Scanner_loc_cell, laser_resolution, beam_divergence, number_scanners)
    %% Rays from each scanner that intersect the tree stem %%
        theta_cell      = cell(1, number_scanners);
        theta_LB_cell   = cell(1, number_scanners);
        theta_UB_cell   = cell(1, number_scanners);
        
        for s = 1 : number_scanners
            % The difference between the scanner and tree
            Scanner_loc         = Scanner_loc_cell{s};
            
            location_difference = Tree_centre - Scanner_loc;
            delta_x             = location_difference(1);
            delta_y             = location_difference(2);
            
            % The angle between the scanner and the tree centre
            viewing_direction   = pi/2 - atan2(delta_y, delta_x);
            
            if viewing_direction < 0
                viewing_direction = 2*pi + viewing_direction;
            end
            
            % Angle alpha at which the beam grazes the tree stem (as seen from the vector between scanner and tree stem) 
            Range = sqrt(sum((location_difference).^2));            % Distance between scanner and tree stem centre
            
            alpha_limit = asin(Tree_radius / Range) + beam_divergence;
            alpha_limit_rounded = laser_resolution * floor(alpha_limit / laser_resolution);

            % The grazing angles, taking the viewing direction into account
            theta_limit_upper   = alpha_limit + viewing_direction;
            theta_limit_lower   = -alpha_limit + viewing_direction;
            
            theta_LB_cell{s}    = theta_limit_lower;
            theta_UB_cell{s}    = theta_limit_upper;
            
            % List of clockwise positive angles as seen from the scanner that hit the tree stem
            theta_list      = (-alpha_limit_rounded : laser_resolution : alpha_limit_rounded)' + viewing_direction;    
            theta_cell{s}   = theta_list;            
        end
end