% This script applies (weighted) Pratt's least squares to fit a circle to the coordinates [x_list, y_list] 
% If w_list is true, they are weighted by [w_list]. Otherwise this input can remain empty

% The step threshold dictates the point at which Newton's method converges, and max iterations the maximum number of iterations that is allowed to take place

% The outputs are the circle centroid [circle_x, circle_y] and its radius circle_r

function [circle_x, circle_y, circle_r] = Pratt_Circle_Fit(x_list, y_list, step_threshold, max_iterations, Weighting, w_list)

    % The number of data points
    n = length(x_list);
    
    % A column vector shape is ensured
    x_list = reshape(x_list, [n, 1]);
    y_list = reshape(y_list, [n, 1]);

    % The weights are equal if weighting is not desired
    if Weighting == false
        w_list = ones(n, 1);
    elseif Weighting == true
        w_list = reshape(w_list, [n, 1]);
    end
        
    % The weighted centroid of the given coordinates
    x_list_w = x_list .* w_list;
    y_list_w = y_list .* w_list;
    
    centroid_x_w = sum(x_list_w) / sum(w_list);
    centroid_y_w = sum(y_list_w) / sum(w_list);

    % The coordinates are centered and weighted
    x_list_c = (x_list - centroid_x_w);
    y_list_c = (y_list - centroid_y_w);

    % The effective radii (z coordinate) at each point is computed
    z_list_c = x_list_c.^2 + y_list_c.^2;

    % The moments the points produce around it in 3D    
    M_xx = sum(x_list_c.^2 .* w_list) / sum(w_list);
    M_xy = sum(x_list_c .* y_list_c .* w_list) / sum(w_list);
    M_xz = sum(x_list_c .* z_list_c .* w_list) / sum(w_list);

    M_yy = sum(y_list_c.^2 .* w_list) / sum(w_list);
    M_yx = sum(y_list_c .* x_list_c .* w_list) / sum(w_list);
    M_yz = sum(y_list_c .* z_list_c .* w_list) / sum(w_list);

    M_zz = sum(z_list_c.^2 .* w_list) / sum(w_list);

    % These are used for the coefficients of the polynomial describing the circle
    M_z = M_xx + M_yy;
    Covariance_xy = M_xx * M_yy - M_xy * M_yx;

    A = M_xz^2*M_yy + M_yz^2*M_xx - M_zz*Covariance_xy - 2*M_xz*M_yz*M_xy + M_z*M_z*Covariance_xy;
    B = M_zz*M_z + 4*Covariance_xy*M_z - M_xz^2 - M_yz^2 - M_z^3;
    C = 4*Covariance_xy - 3*M_z^2 - M_zz;

    % Newton's method
    iter = 0;

    D_new = 0;
    E_new = Inf;

    while iter <= max_iterations
        iter = iter + 1;

        E_old = E_new;
        E_new = A + D_new*(B + D_new*(C + 4*D_new^2));

        % A message if the method fails to converge
        if abs(E_new) > abs(E_old)
            disp('Error: Convergence cannot be obtained')
            D_new = 0;
            break
        end

        F = B + D_new * (C^2 + 16*D_new^2);
        D_old = D_new;
        D_new = D_old - E_new / F;

        % Convergence check
        step = (D_new - D_old) / D_old;

        if abs(step) < step_threshold
            break
        end
    end

    % Check whether the method succeeded
    if iter == max_iterations
        disp('Warning: Newton''s method failed to converge');
    end

    % The circle parameters
    Determinant         = D_new^2 - D_new * M_z + Covariance_xy;

    circle_x_c_w        = (M_xz*(M_yy - D_new) - M_yz*M_xy) / (2*Determinant);
    circle_x            = circle_x_c_w + centroid_x_w;
    circle_y_c_w        = (M_yz*(M_xx - D_new) - M_xz*M_xy) / (2*Determinant);
    circle_y            = circle_y_c_w + centroid_y_w;

    circle_r            = sqrt(sum([circle_x_c_w, circle_y_c_w].*[circle_x_c_w, circle_y_c_w]) + M_z + 2*D_new);
end
