% This script computes the Mahalanobis distance of samples to a Gaussian mixture model, both to the nearest component and using the aggregated multivariate normal
% The distance is expressed in the number of standard deviations from the mean

function [Avg_aggr_Mahal_distance, Avg_min_Mahal_distance] = Mahalanobis_Distance_GMM(data_matrix, aggr_mu_list, aggr_cov_matrix, GM_model, Plot, Print)
    
    %% The Mahalonobis distance %%
        %--% The minimum Mahalonobis distance %--%
        Mahal_distance_matrix   = sqrt(mahal(GM_model, data_matrix));
                        
        % The average minimum distance
        min_Mahal_distance_list = min(Mahal_distance_matrix, [], 2);        % The component to which the data is closest is used
        Avg_min_Mahal_distance  = mean(min_Mahal_distance_list);
        
        %--% The Mahalonobis distance to the weighted aggregate %--%
        number_samples = size(data_matrix, 1);
        
        aggr_Mahal_distance_list = zeros(1, number_samples);
        
        for i = 1 : number_samples
            % The Mahal distance of this sample to the aggregated model
            sample          = data_matrix(i, :);
            Mahal_distance  = sqrt((sample - aggr_mu_list) * (aggr_cov_matrix \ (sample - aggr_mu_list)'));
            
            aggr_Mahal_distance_list(i) = Mahal_distance;
        end
        
        % The average distance to the aggregated model
        Avg_aggr_Mahal_distance = mean(aggr_Mahal_distance_list);
        
    %% Printed result and plot %%    
        if Print == true
            disp('--------------------')
            fprintf('The average Mahalanobis distance to the aggregate model is %g STD \n', Avg_aggr_Mahal_distance);
            fprintf('The average Mahalanobis distance to the nearest component is %g STD \n', Avg_min_Mahal_distance);
            disp('--------------------')
            fprintf('\n')
        end
        
        if Plot == true
            %--% Distances to each component %--%
            % The weight of each component is also displayed
            component_weights_list  = GM_model.ComponentProportion;
            number_components = length(component_weights_list);

            % Colours for the histograms
            hist_colourmap = cbrewer('qual', 'Set2', number_components);
            hist_colourmap = max(hist_colourmap, 0);
            hist_colourmap = min(hist_colourmap, 1);
            
            figure(1)
            
            % Size and white background
            set(gcf, 'Units', 'Normalized', 'Position', [0.05, 0.05, 0.9, 0.85])
            set(gcf, 'color', [1, 1, 1])  
            
            for c = 1 : number_components
                % This subplot's data and colour
                Mahal_distance_list = Mahal_distance_matrix(:, c);
                component_weight    = component_weights_list(c);
                
                hist_colour         = hist_colourmap(c, :);
                
                % The bin-width
                IQR             = iqr(Mahal_distance_list);
                bin_width       = 2 * IQR * length(Mahal_distance_list)^(-1/3);             % Freedman-Diaconis rule
                data_amplitude  = max(Mahal_distance_list) - min(Mahal_distance_list);
                number_bins     = round(data_amplitude / bin_width);
                
                subplot(1, number_components, c);
                hold on
                grid on
                
                % Component number and weight
                title_str = sprintf('Component %g, weight %.2g', c, component_weight);
                title(title_str);
                
                % The histogram
                histogram(Mahal_distance_list, number_bins, 'FaceColor', hist_colour, 'Normalization', 'probability');
                
                % Axes
                xlabel('Mahal. dist. [STD]')
                ylabel('P [-]');
                
                ylim([0, 1]);
                
                % Font size
                set(gca, 'FontSize', 15);
                set(gca, 'LineWidth', 2);
    
                hold off
            end
            
            % The figure is saved and closed
            figure_name = 'Mahalonobis_Distance.png';
            export_fig(1, figure_name);
            
            close(1);
        end
end