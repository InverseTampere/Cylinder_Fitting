% This is an updated version of the script created by Pasi Raumonen for TreeQSM

% ---------------------------------------------------------------------
% Input 
    % cylinder_geometry         Cylinder parameters [x0 y0 alpha beta r]'
    % Geometry_Indices          Structure containing the indices of the geometry parameters to ensure consistency
    % point_cloud_matrix        Point cloud
    % weight_list               (Optional) Weights for the points
     
% Output
    % distance_list             Signed distances of points to the cylinder surface:
    %                           dist(i) = sqrt(xh(i)^2 + yh(i)^2) - r, where 
    %                           [xh yh zh]' = Ry(beta) * Rx(alpha) * ([x y z]' - [x0 y0 0]')
    % Jacobian                  Jacobian matrix d dist(i)/d par(j).
% ---------------------------------------------------------------------

%%%%%% CHANGES %%%%%%
% To ensure consistency, the indices of the geometry variables are required as input
% Changed variable names to be more descriptive
% Changed the Jacobian weighting code to be shorter but slightly less robust (as the weight_list must be given as a colun vector now)
% Simplified the code in a few places
% Corrected the Jacobians for alpha and beta

function [distance_list, sum_squared_distance, Jacobian] = func_grad_cylinder_2(cylinder_geometry, Geometry_Indices, point_cloud_matrix, weight_list)

    %% Cylinder parameters %%
        % Five cylinder parameters: 
        % Location = axis point intersects xy-plane: x0 and y0 (z0 == 0)
        % Rotation angles around x and y axis = alpha and beta
        % Radius = r

        [centre_ind, alpha_ind, beta_ind, radius_ind]   = deal(Geometry_Indices.centre, Geometry_Indices.alpha, Geometry_Indices.beta, Geometry_Indices.radius);
        [proj_centre, alpha, beta, radius]              = deal(cylinder_geometry(centre_ind), cylinder_geometry(alpha_ind), cylinder_geometry(beta_ind), cylinder_geometry(radius_ind));
    
        % The projected centre is made three dimensional
        proj_centre_3D = [proj_centre; 0]';

    %% Point cloud transformation %%
        % Transformed points:
        % Pt = [xh yx zh] = Ry(beta) * Rx(alpha) * (P - [x0 y0 0])
    
        % The rotation matrices and their derivatives
        [R, R_alpha, R_beta, DR_alpha, DR_beta] = form_rotation_matrices_2(alpha, beta);
    
        % Translation and rotation
        point_cloud_matrix_t = (point_cloud_matrix - proj_centre_3D) * R';

        % "Plane points":
        % Qt = Pt * I2 = [xh yh];
        [x_t_list, y_t_list] = deal(point_cloud_matrix_t(:, 1), point_cloud_matrix_t(:, 2));

    %% Distances %%
        % Distance:
        % D(x0,y0,alpha,beta,r) = sqrt( dot(Qt,Qt) )-r = sqrt( Qt*Qt' )-r
        % rt = sqrt( dot(Qt,Qt) )

        % Calculate the distances
        radius_list     = sqrt(x_t_list.^2 + y_t_list.^2);
        distance_list   = radius_list - radius;                 % Distances to the cylinder surface
        
        % Weighted, if the weights are given
        if nargin == 4
            distance_list = weight_list .* distance_list;      
        end

        % Sum of squared distances
        sum_squared_distance = norm(distance_list);

    %% Jacobian %%
        if nargout > 1
            % Initialise the Jacobian
            num_points          = size(point_cloud_matrix, 1);
            num_geom_parameters = length(cylinder_geometry);

            Jacobian = zeros(num_points, num_geom_parameters);
            
            % Normalised coordinates (N = Qt / rt)
            N = [x_t_list./radius_list, y_t_list./radius_list];

            % Derivatives w.r.t. the projected centre
            neg_x_axis      = [-1, 0, 0];
            rot_neg_x_axis  = (R * neg_x_axis')';

            Jacobian(:, centre_ind(1)) = sum(rot_neg_x_axis(1:2) .* N, 2);

            neg_y_axis      = [0, -1, 0];
            rot_neg_y_axis  = (R * neg_y_axis')';

            Jacobian(:, centre_ind(2)) = sum(rot_neg_y_axis(1:2) .* N, 2);
            
            % Derivative w.r.t. alpha (rotation around x-axis)
            A = point_cloud_matrix_t * DR_alpha' * R_beta';
            
            Jacobian(:, alpha_ind) = sum(A(:, 1:2) .* N, 2);
            
            % Derivative w.r.t. beta (rotation around y-axis)
            B = point_cloud_matrix_t * R_alpha' * DR_beta';
            
            Jacobian(:, beta_ind) = sum(B(:, 1:2) .* N, 2);
            
            % Derivative w.r.t. the radius is simply -1
            Jacobian(:, radius_ind) = repmat(-1, [num_points, 1]);

            % Weighted Jacobian
            if nargin == 4
                Jacobian = weight_list .* Jacobian;
            end            
        end    
end
