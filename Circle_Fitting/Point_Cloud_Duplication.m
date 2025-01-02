% This script duplicates the point cloud on the other side of the circle, to avoid bias

function [x_duplicated_points_cell, y_duplicated_points_cell] = Point_Cloud_Duplication(Input_Geometry, x_points_cell, y_points_cell, number_scanners)

    %% The points of each scanner are duplicated with respect to the circle centre %%
        circle_x                    = Input_Geometry(1);
        circle_y                    = Input_Geometry(2);
    
        x_duplicated_points_cell    = cell(1, number_scanners);
        y_duplicated_points_cell    = cell(1, number_scanners);
    
        for s = 1 : number_scanners
            % This scanner's point cloud
            x_points_list = x_points_cell{s};
            y_points_list = y_points_cell{s};
                        
            % Additionally, points are duplicated on the other side, to avoid bias            
            centroid_vector_x       = circle_x - x_points_list;                 % Vector toward the centroid of the circle
            centroid_vector_y       = circle_y - y_points_list;
            
            x_mirror_points_list    = x_points_list + 2 * centroid_vector_x;    % The mirrored points are located twice this vector away
            y_mirror_points_list    = y_points_list + 2 * centroid_vector_y;
            
            x_total_points_cell     = {x_points_list, x_mirror_points_list};
            y_total_points_cell     = {y_points_list, y_mirror_points_list};
            
            x_duplicated_points_list    = vertcat(x_total_points_cell{:});      % Doesn't serve a purpose other than to avoid the unecessary warning
            y_duplicated_points_list    = vertcat(y_total_points_cell{:});
            
            x_duplicated_points_cell{s} = x_duplicated_points_list;
            y_duplicated_points_cell{s} = y_duplicated_points_list;
        end
end