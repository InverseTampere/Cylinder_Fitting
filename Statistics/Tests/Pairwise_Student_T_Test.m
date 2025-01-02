% The pairwise Student's t-test is performed for the variables in the given structure arrays. P-values are returned in a structure containing the same variables

function P_Values = Pairwise_Student_T_Test(Structure_A, Structure_B)
    
    %% P-values resulting from Student's t-test %%
        % Assessed variables
        variable_labels     = fieldnames(Structure_A(1));
        number_variables    = length(variable_labels);
    
        % Student's t-test
        P_Values = struct();

        for v = 1 : number_variables
            % The structures' data
            variable_label      = variable_labels{v};
            
            structure_A_cell    = {Structure_A.(variable_label)};
            structure_A_matrix  = vertcat(structure_A_cell{:});
            structure_B_cell    = {Structure_B.(variable_label)};
            structure_B_matrix  = vertcat(structure_B_cell{:});

            % If the data is multi-dimensional, each dimension is assessed indepently
            num_dim         = size(structure_A_matrix, 2);
            variable_P_list = zeros(1, num_dim);

            for d = 1 : num_dim
                % This dimension's data
                structure_A_list = structure_A_matrix(:, d);
                structure_B_list = structure_B_matrix(:, d);

                % Student's t-test
                [~, P]              = ttest(structure_A_list, structure_B_list);
                variable_P_list(d)  = P;
            end

            % The P-values are added to the structure
            P_Values.(variable_label) = variable_P_list;
        end
end