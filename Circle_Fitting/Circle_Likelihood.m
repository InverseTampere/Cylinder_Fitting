% This function computes the normalised likelihood as a function of the angle around the circle, theta
% It revolves around a translation to the [r, p] coordinate frame, where r is in the radial direction and p in the propagation direction

function Norm_Likelihood_list = Circle_Likelihood(theta_list, x_point, y_point, scanner_x, scanner_y, circle_x, circle_y, circle_radius, sigma_rad, sigma_prop, laser_radius)
    %% Translation of the circle centre to the r,p coordinate frame %%    
        % The angle from the scanner to the point
        delta_x = x_point - scanner_x;
        delta_y = y_point - scanner_y;

        alpha   = pi/2 - atan2(delta_y, delta_x);

        circle_centre_star = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)] * [circle_x - x_point; circle_y - y_point];
        circle_r = circle_centre_star(1);
        circle_p = circle_centre_star(2);
        
    %% The angles covered by the beamwidth %%
        % Angles that are not covered have a likelihood of zero. They are now given the NaN flag
        theta_a = asin((sign(circle_r) * laser_radius - circle_r) / circle_radius);
        
        % In case the angle is negative, it is made positive again
        if theta_a < 0
            theta_a = pi - theta_a;
        end

        if -sign(circle_r)*circle_r + circle_radius < laser_radius    % The laser beam lies partially outside the circle
            theta_b = (2 + sign(circle_r))*pi - theta_a;
            
            theta_list(theta_list < theta_a) = NaN;
            theta_list(theta_list > theta_b) = NaN;
        else                                                % Otherwise there are two separate areas where the laser intersects the circle
            theta_b = -asin((sign(circle_r) * laser_radius + circle_r) / circle_radius);
            
            % In case the angle is negative, it is made positive again
            if theta_b < 0
                theta_b = pi - theta_b;
            end
            
            theta_d = (2 + sign(circle_r))*pi - theta_a;
            theta_c = theta_d - (theta_b - theta_a);

            theta_list(theta_list < theta_a) = NaN;
            theta_list(theta_list > theta_b & theta_list < theta_c) = NaN;
            theta_list(theta_list > theta_d) = NaN;            
        end
            
    %% The normalised likelihood function %%
        G = -((circle_r/sigma_rad)^2 + (circle_p/sigma_prop)^2)/2;
        a = -(circle_radius / sigma_rad)^2 / 2;
        b = -(circle_radius / sigma_prop)^2 / 2;
        c = -circle_radius * circle_r / sigma_rad^2;
        d = -circle_radius * circle_p / sigma_prop^2;

        Norm_Likelihood_Fun = @(theta) exp(a*sin(theta).^2 + b*cos(theta).^2 + c*sin(theta) + d*cos(theta) + G);
    
    %% The normalised likelihood %%
        Norm_Likelihood_list = Norm_Likelihood_Fun(theta_list);
        
        % The likelihoods outside of the laser beam have a NaN value, which is substituted for zero
        Norm_Likelihood_list(isnan(Norm_Likelihood_list)) = 0;
end