% The point cloud, distributions and scanner locations may be translated w.r.t. the given point
% Undesired inputs may be left empty

function [Point_Cloud_Coord_t, Point_Cloud_Distributions_t, point_cloud_matrix_t, point_cloud_centroid_t, Scanning_Parameters_t, Scanner_loc_cell_t, number_scanners] = Point_Cloud_Translation(translation_point, Point_Cloud_Coord, Point_Cloud_Distributions, Scanning_Parameters)

    %% Point cloud centering %% 
        if ~isempty(Point_Cloud_Coord)
            % The point cloud
            x_points_cell       = Point_Cloud_Coord.x_points_cell;
            y_points_cell       = Point_Cloud_Coord.y_points_cell;
            z_points_cell       = Point_Cloud_Coord.z_points_cell;

            point_cloud_matrix  = [vertcat(x_points_cell{:}), vertcat(y_points_cell{:}), vertcat(z_points_cell{:})];
    
            % Centering the point cloud
            point_cloud_matrix_t    = point_cloud_matrix - translation_point;
            point_cloud_centroid_t  = mean(point_cloud_matrix_t, 1);

            Point_Cloud_Coord_t                         = Point_Cloud_Coord;
            Point_Cloud_Coord_t.Point_cloud_matrix      = point_cloud_matrix_t;
            Point_Cloud_Coord_t.Point_cloud_centroid    = point_cloud_centroid_t;
            
            Translation_x_fun                   = @(x) x - translation_point(1);
            x_points_cell_t                     = cellfun(Translation_x_fun, x_points_cell, 'UniformOutput', false);
            Point_Cloud_Coord_t.x_points_cell   = x_points_cell_t;
    
            Translation_y_fun                   = @(y) y - translation_point(2);
            y_points_cell_t                     = cellfun(Translation_y_fun, y_points_cell, 'UniformOutput', false);
            Point_Cloud_Coord_t.y_points_cell   = y_points_cell_t;
    
            Translation_z_fun                   = @(z) z - translation_point(3);
            z_points_cell_t                     = cellfun(Translation_z_fun, z_points_cell, 'UniformOutput', false);
            Point_Cloud_Coord_t.z_points_cell   = z_points_cell_t;        
        
        else
            % An empty coordinate output is returned
            Point_Cloud_Coord_t     = [];
        end

    %% Distribution centering %%
        if ~isempty(Point_Cloud_Distributions)
            % The expected value of each distribution is used as the point cloud
            distribution_mu_cell    = Point_Cloud_Distributions.distribution_mu_cell;
            point_cloud_matrix      = vertcat(distribution_mu_cell{:});

            point_cloud_matrix_t    = point_cloud_matrix - translation_point;
            point_cloud_centroid_t  = mean(point_cloud_matrix_t, 1);
            
            % The expected values are shifted by the point
            Translation_fun         = @(distribution_mu) distribution_mu - translation_point;
            distribution_mu_cell_t  = cellfun(Translation_fun, distribution_mu_cell, 'UniformOutput', false);

            % The new structure
            Point_Cloud_Distributions_t                         = Point_Cloud_Distributions;
            Point_Cloud_Distributions_t.distribution_mu_cell    = distribution_mu_cell_t;
        else
            Point_Cloud_Distributions_t = [];
        end

    %% Centering the scanner locations %%
        if ~isempty(Scanning_Parameters)
            % The scanner locations
            Scanner_loc_cell    = Scanning_Parameters.Scanner_loc_cell;
            number_scanners     = Scanning_Parameters.number_scanners;
            
            % Translation
            Scanner_Translation_fun = @(scanner_loc) scanner_loc - translation_point;
            Scanner_loc_cell_t      = cellfun(Scanner_Translation_fun, Scanner_loc_cell, 'UniformOutput', false);

            Scanning_Parameters_t   = struct('Scanner_loc_cell', {Scanner_loc_cell_t}, 'number_scanners', number_scanners);
        else
            [Scanning_Parameters_t, Scanner_loc_cell_t, number_scanners] = deal([]);
        end
end