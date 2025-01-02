% This script makes the given matrix as square as possible through permutation
% Symmetry is quantified as the sum of the matrix minus its transpose
% If the given matrix is not square, it is made square by adding zero-entries

function [symmetric_matrix, symmetric_sum] = Matrix_Symmetry(matrix)

    %% Ensuring squareness %%
        % Size of the matrix
        [number_rows, number_columns] = size(matrix);

        % Adding zero-entries to make the matrix square
        matrix = [matrix, zeros(number_rows, number_rows - number_columns)];
        matrix = [matrix; zeros(number_columns - number_rows, number_columns)];
        
        [~, number_columns] = size(matrix);

    %% Assessing symmetry %%
        % Possible permutations of the columns of the matrix
        perm_matrix     = perms(1 : number_columns);
        number_perms    = factorial(number_columns);
        perm_cell       = mat2cell(perm_matrix, ones(1, number_perms), number_columns);

        % The permuted matrices
        Permutation_fun     = @(column_order) matrix(:, column_order);
        perm_matrix_cell    = cellfun(Permutation_fun, perm_cell, 'UniformOutput', false);

        % Their symmetry measure
        Symmetry_Sum_fun    = @(matrix) sum(matrix - matrix', 'all');
        symmetry_sum_list   = cellfun(Symmetry_Sum_fun, perm_matrix_cell);

        % The most symmetric matrix
        [symmetric_sum, symm_ind]   = min(symmetry_sum_list);
        symmetric_matrix            = perm_matrix_cell{symm_ind};

end