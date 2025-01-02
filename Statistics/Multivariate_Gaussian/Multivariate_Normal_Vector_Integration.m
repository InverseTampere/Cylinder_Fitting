% This script determines the vector integral of a 3D multivariate normal probability density function using geometric projection
% The vector integral is taken along the vector
% The normalised likelihood values are the density values, but made independent of dimension and uncertainty

function [f_vector_list, L_N_vector_list] = Multivariate_Normal_Vector_Integration(proj_points_matrix, alpha, proj_ellipsoid_sigmae, Plot)

    %% Coverage probability %%
        % The number of standard deviations follows from the coverage probability
        num_dim = size(proj_points_matrix, 2);
        P       = 1 - alpha;
        m_STD   = sqrt(chi2inv(P, num_dim));

    %% Integration onto the vector %%
        % Height w.r.t. the projection
        height_list     = proj_points_matrix(:, num_dim);
        height_sigma    = proj_ellipsoid_sigmae(num_dim);
        
        % Normalised likelihood values
        L_N_vector_list         = exp(-1/2 * (height_list / height_sigma).^2);
        
        % Probability density values
        normalisation_factor    = 1/(sqrt(2*pi) * height_sigma);
        dimension_factor        = P / chi2cdf(m_STD^2, 1);              % As the Gaussian is now considered 1-dimensional
        f_vector_list           = normalisation_factor * dimension_factor * L_N_vector_list;
        
    %% Plot %%
        if Plot == true
            % Message showing the integrated density value
            delta_height = max(height_list) - min(height_list);
            fprintf('The integrated density over the vector: %.3g \n', mean(f_vector_list) * delta_height);
            
            %--% Density plot %--%
            figure(1)
            % Set the size and white background color
            set(gcf, 'Units', 'Normalized', 'Position', [0.05 0.05 0.9 0.85])
            set(gcf, 'color', [1, 1, 1])     
  
            hold on
            grid on
            
            f_UB = 1.2 * max([f_vector_list; 1e-6]);
            
            plot(height_list, f_vector_list, 'color', 'r', 'LineWidth', 2, 'DisplayName', 'f(h)')
            plot([0, 0], [0, f_UB], 'color', 'k', 'LineWidth', 2, 'DisplayName', '\mu_{h}');
            
            % Axes
            xlabel('h [m]');
            ylabel('f (m^{-2})');
            ylim([0, f_UB]);
            
            % Legend
            legend('show', 'location', 'northoutside');

            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);
            
            % Pause
            disp('The script will continue and the figure will close upon a button-press');
            pause();
            close(1);                
        end
        
end