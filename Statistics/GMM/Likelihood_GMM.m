% This script computes the likelihood of the Gaussian mixture model, given the input data

function [L_norm_avg, L_norm_max, L_norm_min, L_norm_list] = Likelihood_GMM(Data_matrix, GM_model, Print, Plot)
    
    %% The Gaussian mixture model properties %%
        number_components   = GM_model.NumComponents;           % The number of multivariate normal distributions (components)
        mean_vector_matrix  = GM_model.mu;                      % The mean vectors of the components
        Covariance_matrix   = GM_model.Sigma;                   % The covariance matrices
        Shared_Covariance   = GM_model.SharedCovariance;        % Whether or not each component has its own covariance matrix
        weights_list        = GM_model.ComponentProportion;     % The weight of each component (sums to 1)
        
    %% The normalised likelihood %%
        [num_samples, num_dim] = size(Data_matrix);
        
        % The normalised likelihood is computed for all samples per component
        L_norm_component_matrix = zeros(num_samples, number_components);
        
        for k = 1 : number_components
            % This component's properties
            mu_list = mean_vector_matrix(k, :);
            
            if Shared_Covariance == true
                Sigma = Covariance_matrix;
            else
                Sigma = Covariance_matrix(:, :, k);
            end
            
            % The normalised likelihood values of each sample
            for i = 1 : num_samples
                sample = Data_matrix(i, :);
                
                % The Mahalanobis distance between the sample and this component
                Mahal_distance  = sqrt((sample - mu_list) * (Sigma \ (sample - mu_list)'));
                
                % The normalised likelihood value, corrected for dimension
                L_norm = sqrt(2)^num_dim * exp(-1/2 * Mahal_distance^2);
                
                L_norm_component_matrix(i, k) = L_norm;
            end
        end
        
        % The values are weighted and summed over all components
        L_norm_component_matrix_w   = weights_list .* L_norm_component_matrix;
        L_norm_list                 = sum(L_norm_component_matrix_w, 2);
        
        % The average and extrema
        L_norm_avg = mean(L_norm_list);
        L_norm_max = max(L_norm_list);
        L_norm_min = min(L_norm_list);
        
    %% Printed result and histogram %%
        if Print == true
            disp('--------------------')
            fprintf('The average normalised likelihood of the GMM: %.3g \n', L_norm_avg); 
        end
        
        if Plot == true
            figure(1)
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.85])
            set(gcf, 'color', [1, 1, 1])     

            hold on
            grid on

            % Histogram
            IQR             = iqr(L_norm_list);
            bin_width       = 2 * IQR * length(L_norm_list)^(-1/3);     % Freedman-Diaconis rule
            data_amplitude  = max(L_norm_list) - min(L_norm_list);
            number_bins     = round(data_amplitude / bin_width);
            histogram(L_norm_list, number_bins, 'FaceColor', 'r', 'normalization', 'probability');

            % Axes
            xlabel('avg. L_N [-]');
            ylabel('P [-]');

            ylim([0, 1]);
            xlim([0, 1]);

            % Font size
            set(gca, 'FontSize', 15);
            set(gca, 'LineWidth', 2);

            hold off
            
            % Figure is saved and closed
            export_fig(1, 'GMM_Likelihood_Histogram.png');
            
            close(1);
        end
    
end