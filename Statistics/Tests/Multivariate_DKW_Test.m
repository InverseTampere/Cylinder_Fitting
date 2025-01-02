% This script applies the Dvoretzky-Kiefer-Wolfowitz inequality to compute the bound around a multivariate CDF (Naaman, 2021 [A108])
% The hypothesis that the data fits the Gaussian mixture model with significance alpha is marked by H

function [H, alpha_value, KS_value] = Multivariate_DKW_Test(data_matrix, GM_model, alpha, Print)
    %% The Kolmogorov-Smirnov statistic value %%
    
        % As the CDF is monotonic, the supremum must lie at a sample location or that of a line intersect (Gosset, 1987 [A107])

        %--% Samples %--%
        [num_samples, num_dim] = size(data_matrix);
        
        D_samples_list = zeros(1, num_samples);
        
        for i = 1 : num_samples
            % The sample data
            sample = data_matrix(i, :);
            
            % The true CDF value
            F = cdf(GM_model, sample);
            
            % The empirical CDF value
            sample_diff_matrix  = data_matrix < sample;
            sample_diff_list    = sum(sample_diff_matrix, 2);

            num_smaller_samples = length(find(sample_diff_list == num_dim));    % Smaller in all dimensions  
        
            EF_greater  = (num_smaller_samples + 1) / num_samples;              % A point slightly greater than this sample location, i.e. it includes the current sample
            EF_smaller  = num_smaller_samples / num_samples;                    % It does not include the current sample, as it is slightly smaller
            
            % The largest absolute difference value is used
            D_greater   = abs(EF_greater - F);
            D_smaller   = abs(EF_smaller - F);
            
            D_samples_list(i) = max(D_greater, D_smaller);
        end
        
        %--% Line intersects %--%
        D_intersects_list = zeros(1, num_samples);
        
        for i = 1 : num_samples
            % The sample data
            sample = data_matrix(i, :);
            
            % Points with which intersection can take place are biger in one dimension, and smaller in all others
            offset_sign_matrix  = data_matrix < sample;
            offset_sign_list    = sum(offset_sign_matrix, 2);
            candidate_points    = offset_sign_list == num_dim - 2;      % 2 is the difference between positive and negative sign
            candidate_points(i) = false;                                % The current point has a sign offset of 0, which can falsely be seen as a candidate if num_dim is 2

            % The intersects are then the maximum in each dimension
            data_candidate_points   = data_matrix(candidate_points, :);
            data_intersects         = max(data_candidate_points, sample);
            
            % The absolute difference of each intersect
            number_intersects           = size(data_intersects, 1);
            
            D_sample_intersects_list    = zeros(1, number_intersects);
            
            for j = 1 : number_intersects
                % This intersect's data
                intersect = data_intersects(j, :);
                
                % The true CDF value
                F = cdf(GM_model, intersect);
                
                % The empirical CDF value
                sample_diff_matrix  = data_matrix < intersect;
                sample_diff_list    = sum(sample_diff_matrix, 2);

                num_smaller_samples = length(find(sample_diff_list == num_dim));            
        
                EF_greater  = (num_smaller_samples + 1) / num_samples;       % A point slightly greater than this sample location, i.e. it includes the current sample
                EF_smaller  = num_smaller_samples / num_samples;             % It does not include the current sample, as it is slightly smaller
            
                % The largest absolute difference value is used
                D_greater   = abs(EF_greater - F);
                D_smaller   = abs(EF_smaller - F);

                D_sample_intersects_list(j) = max(D_greater, D_smaller);
            end
            
            % The maximum of the intersects is appended
            if number_intersects > 0
                D_intersects_list(i) = max(D_sample_intersects_list);
            end
        end
        
        %--% The Kolmogorov-Smirnov statistic value %--%
        D_total_list    = [D_intersects_list, D_samples_list];
        KS_value        = max(D_total_list);
        
    %% The statistical significance %%
        % The probability that the supremum is greater than the found value, is bounded by the multivariate DKW inequality
        alpha_value = 2 * num_dim * exp(-2 * num_samples * KS_value^2);         % Note that the number of samples is assumed to be sufficiently large
        
        % The test is rejected if it is greater then the given critical value        
        if alpha_value > alpha
            H = false;
        else
            H = true;
        end
        
        if Print == true
            disp('--------------------')
            if H == false
                fprintf('The DKW inequality is NOT met, alpha = %.3g > alpha_cr = %.3g \n', alpha_value, alpha);
            else
                fprintf('The DKW inequality is met, alpha = %.3g < alpha_cr = %.3g \n', alpha_value, alpha);
            end
        end
        
end