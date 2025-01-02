% Translation of the coordinates from a point on the circle to the radial and propagation offsets based on the hit and the direction of the beam

function [r, p] = Coordinate_Translator(theta, circle_radius, circle_centroid_x, circle_centroid_y, x_point, y_point, scanner_loc)
    % Coordinates on the circle, w.r.t. the circle centroid
    xc = circle_radius * sin(theta);
    yc = circle_radius * cos(theta);
    
    % Translation to the location of the point
    xc_t = xc - x_point + circle_centroid_x;
    yc_t = yc - y_point + circle_centroid_y;
    
    % Angle from the vertical to the laser beam
    delta_x = x_point - scanner_loc(1);
    delta_y = y_point - scanner_loc(2);

    alpha   = pi/2 - atan2(delta_y, delta_x);    
            
    % Rotation to the (r, p) coordinate frame
    r = cos(alpha) * xc_t - sin(alpha) * yc_t;
    p = sin(alpha) * xc_t + cos(alpha) * yc_t;    
end