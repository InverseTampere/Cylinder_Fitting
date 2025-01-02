% This script computes the Bhattacharyya distance and coefficient between two Gaussian mixture models (GMMs)

% The weighted sum of the pair-wise distances according to Lahlimi et al. 2017 [A119] is nonzero even for identical GMMs, 
% so it is instead computed as the sum of the matrix minus its transpose which acts as a symmetry measure

function [Bhatt_coeff, Bhatt_distance] = Bhattacharyya_GMM(GM_Model_a, GM_Model_b)

    %% GMM properties %%
        % Model a (k is the number of Gaussians, m the number of variables)
        number_components_a         = GM_Model_a.NumComponents;             % Number of Gaussians in the model
        Mu_matrix_a                 = GM_Model_a.mu;                        % Expected values (k x m)
        Sigma_matrix_a              = GM_Model_a.Sigma;                     % Covariance matrix (m x m x k or if covariance is shared m x m)
        Component_weights_list_a    = GM_Model_a.ComponentProportion;       % The weight of each Gaussian (1 x k)
        Shared_Covariance_a         = GM_Model_a.SharedCovariance;          % Whether or not the components have the same or different covariance matrices

        % Model b (k is the number of Gaussians, m the number of variables)
        number_components_b         = GM_Model_b.NumComponents;             % Number of Gaussians in the model (k)
        Mu_matrix_b                 = GM_Model_b.mu;                        % Expected values (k x m)
        Sigma_matrix_b              = GM_Model_b.Sigma;                     % Covariance matrix (m x m x k or if covariance is shared m x m)
        Component_weights_list_b    = GM_Model_b.ComponentProportion;       % The weight of each Gaussian (1 x k)
        Shared_Covariance_b         = GM_Model_b.SharedCovariance;          % Whether or not the components have the same or different covariance matrices

    %% Pair-wise Bhattacharyya distances %%
        % The Bhattacharyya distance between each component of models a and b
        B_distance_matrix = zeros(number_components_a, number_components_b);

        for i = 1 : number_components_a
            % The i'th component's properties of model a
            mu_i        = Mu_matrix_a(i, :);

            if Shared_Covariance_a == true
                Sigma_i = Sigma_matrix_a;
            else
                Sigma_i = squeeze(Sigma_matrix_a(:, :, i));
            end

            for j = 1 : number_components_b
                % The j'th component's properties of model b
                mu_j        = Mu_matrix_b(j, :);
    
                if Shared_Covariance_b == true
                    Sigma_j = Sigma_matrix_b;
                else
                    Sigma_j = squeeze(Sigma_matrix_b(:, :, j));
                end

                % Their bhattacharyya distance
                Sigma       = (Sigma_i + Sigma_j) / 2;
                B_distance  = 1/8*((mu_i - mu_j) / Sigma) * (mu_i - mu_j)' + 1/2*log(det(Sigma)/sqrt(det(Sigma_i)*det(Sigma_j)));            
                
                B_distance_matrix(i, j) = B_distance;
            end
        end

        % The distances are weighted by the component weights
        B_distance_matrix_weighted = Component_weights_list_a' .* Component_weights_list_b .* B_distance_matrix;

    %% Over-all Bhattacharyya distance and coefficient %%
        % The Bhattacharyya distance of the most symmetrical matrix minus its transpose
        [~, Bhatt_distance] = Matrix_Symmetry(B_distance_matrix_weighted);

        % The Bhattacharyya coefficient
        Bhatt_coeff = exp(-Bhatt_distance);

end