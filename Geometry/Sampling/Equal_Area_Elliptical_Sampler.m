% This script creates samples that are equally spaced and covering an elliptical area
% The total number of samples is lower than the initially given number of samples, as an exact ellipse discretisation requires an infinite number of samples

function [ellipse_samples_matrix, number_samples, delta_area] = Equal_Area_Elliptical_Sampler(ellipse_centre, ellipse_radii, ellipse_axes, number_samples)
        
    %% Samples of the ellipse %%
        % The number of samples along each axis for equal distance sampling
        number_alpha_samples    = sqrt(4*number_samples/pi * ellipse_radii(1)/ellipse_radii(2));
        number_alpha_samples    = ceil(number_alpha_samples);
        number_beta_samples     = sqrt(4*number_samples/pi * ellipse_radii(2)/ellipse_radii(1));
        number_beta_samples     = ceil(number_beta_samples);
        
        alpha_steps_list        = ellipse_radii(1) * linspace(-1, 1, number_alpha_samples);
        beta_steps_list         = ellipse_radii(2) * linspace(-1, 1, number_beta_samples);
        
        [alpha_steps_matrix, beta_steps_matrix]  = meshgrid(alpha_steps_list, beta_steps_list); 
        
        % The steps that are inside the ellipse
        ellipse_ind     = (alpha_steps_matrix/ellipse_radii(1)).^2 + (beta_steps_matrix/ellipse_radii(2)).^2 <= 1;
        number_samples  = sum(sum(ellipse_ind));
        
        alpha_list  = alpha_steps_matrix(ellipse_ind);
        beta_list   = beta_steps_matrix(ellipse_ind);
        
        % The resulting samples, rotated by the ellipse axes       
        alpha_axis  = ellipse_axes(1, :);
        beta_axis   = ellipse_axes(2, :);
        
        ellipse_samples_matrix = ellipse_centre + alpha_list .* alpha_axis + beta_list .* beta_axis;     

        % The area per sample
        delta_area = pi * prod(ellipse_radii) / number_samples;
        
end