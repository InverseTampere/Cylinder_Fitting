% This script fits the optimal Gaussian mixture model (GMM) to the given variable data
% If one could not be fitted, an empty model is given. Other outputs are NaN

function [Gaussian_Mixture_Model, Shared_Covariance, number_GM_components, AICc_min] = Gaussian_Mixture_Model_Fitting(variable_matrix, number_GMM_fitment_iter)

    %% Gaussian mixture model fitting %%
        % Fitment options
        GMM_Options = statset('MaxIter', number_GMM_fitment_iter);

        % The maximum number of components may be limited by the numer of samples
        [number_samples, number_variables] = size(variable_matrix);

        max_components = floor(number_samples / (number_variables * (number_variables + 1) / 2 + number_variables + 1));
        max_components = min(max_components, 3 * number_variables);        % More than three times the number of variables is not desired

        number_GM_components_list   = 1 : max_components;
        number_tests                = length(number_GM_components_list);

        %--% Fitment tests %--%
        % For each test, the AIC and model are saved
        AICc_list = NaN(1, number_tests);
        GMM_Tests = struct('GMM', []);

        for t = 1 : number_tests
            number_GM_components = number_GM_components_list(t);

            try
                % In case the covariance matrix is ill-conditioned, it is made equal for every component 
                try
                    Shared_Covariance = false;
                    GMM = fitgmdist(variable_matrix, number_GM_components, 'SharedCovariance', Shared_Covariance, 'Options', GMM_Options);
                catch
                    Shared_Covariance = true;
                    GMM = fitgmdist(variable_matrix, number_GM_components, 'SharedCovariance', Shared_Covariance, 'Options', GMM_Options);
                end

                % It is possible that due to rounding errors a covariance matrix is not positive semi-definite (i.e. the determinant is negative)
                Sigma_matrix = GMM.Sigma;
                    
                if number_GM_components > 1 && Shared_Covariance == false
                    Sigma_det_list  = zeros(1, number_GM_components);
    
                    for g = 1 : number_GM_components
                        Sigma               = squeeze(Sigma_matrix(:, :, g));
                        Sigma_det_list(g)   = det(Sigma);
                    end
                else
                    Sigma_det_list = det(Sigma_matrix);
                end

                if min(Sigma_det_list) < 0              % If it is negative, the script continues with a NaN for this model's AICc
                    continue
                end

                GMM_Tests(t).GMM = GMM;

                % The AIC, corrected for sample size
                AICc            = GMM_Corrected_AIC(GMM, Shared_Covariance, number_samples, number_variables, number_GM_components);
                AICc_list(t)    = AICc;

                % If the AIC has increased, continuing is unecessary
                if t > 1
                    if AICc_list(t) > AICc_list(t - 1)
                        break
                    end
                end

            catch
                % If this also fails, the script continues and the AIC remains NaN
            end
        end

        %--% The optimal number of components is selected using the AIC %--%
        % The model is used which has the minimal (valid) AICc value
        [AICc_min, Opt_model_ind]   = min(AICc_list);

        if ~isnan(AICc_min)  
            Gaussian_Mixture_Model  = GMM_Tests(Opt_model_ind).GMM;
            number_GM_components    = number_GM_components_list(Opt_model_ind);
            Shared_Covariance       = Gaussian_Mixture_Model.SharedCovariance;
            
        % Otherwise empty/NaN outputs are given
        else
            Gaussian_Mixture_Model = [];
            [Shared_Covariance, number_GM_components, AICc_min] = deal(NaN);
        end

    %% Local function to compute the corrected AIC %%
        function AICc = GMM_Corrected_AIC(Gaussian_Mixture_Model, Shared_Covariance, number_samples, number_variables, number_GM_components) 
          % The number of parameters
            mu_par      = number_variables * number_GM_components;
            weight_par  = number_variables;

            if Shared_Covariance == true        % The covariance matrix is shared
                sigma_par   = number_variables * (number_variables + 1) / 2;    % Note that the covariance matrix is diagonally symmetric
            else
                sigma_par   = number_variables * (number_variables + 1) / 2 * number_GM_components;
            end

            num_par     = mu_par + sigma_par + weight_par;

            % The AIC is corrected to account for sample size
            AIC         = Gaussian_Mixture_Model.AIC;            
            AIC_penalty = (2*num_par^2 + 2*num_par) / (number_samples - num_par - 1);

            if AIC_penalty > 0
                AICc    = AIC + AIC_penalty;
            else
                AICc    = NaN;
            end
        end
end