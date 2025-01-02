% This script computes the gradient of the given matrix for all its dimensions
% Forward difference:   FD = (F(x_(i + 1)) - F(x_i)) / dx           Last entry is NaN
% Backward difference:  BD = (F(x_i) - F(x_(i - 1))) / dx           First entry is NaN
% Central difference:   CD = (FD + BD) / 2                          First and last entries are NaN

function [Forward_difference_cell, Backward_difference_cell, Central_difference_cell] = Matrix_Gradient(Matrix, step_size_list)

    %% The gradient of each dimension %%
        % The size of the matrix and the number of dimensions
        matrix_size_list    = size(Matrix);
        number_dimensions   = length(matrix_size_list);
        
        % Matrices of the gradients are saved for each dimension in a cell array
        Forward_difference_cell     = cell(1, number_dimensions);
        Backward_difference_cell    = cell(1, number_dimensions);
        Central_difference_cell     = cell(1, number_dimensions);
        
        for d = 1 : number_dimensions
            % This dimension's (step) size
            dimension_size  = matrix_size_list(d);
            step_size       = step_size_list(d);

            if dimension_size == 1      % A gradient cannot be computed
                forward_difference_matrix   = NaN(matrix_size_list);
                backward_difference_matrix  = NaN(matrix_size_list);
                central_difference_matrix   = NaN(matrix_size_list);

            else
                % The indices for the current and next step
                indices_current     = 1 : dimension_size - 1;
                indices_next        = 2 : dimension_size;

                % Their data
                data_next           = Matrix_Data_Selector(Matrix, d, number_dimensions, indices_next);
                data_current        = Matrix_Data_Selector(Matrix, d, number_dimensions, indices_current);

                % Entries that cannot be computed are given a NaN value
                NaN_size_list       = matrix_size_list;
                NaN_size_list(d)    = 1;
                NaN_entries         = NaN(NaN_size_list);

                % Forward, backward and central differences
                data_delta_list     = (data_next - data_current) / step_size;
                
                forward_difference_matrix   = cat(d, data_delta_list, NaN_entries);
                backward_difference_matrix  = cat(d, NaN_entries, data_delta_list);
                central_difference_matrix   = (forward_difference_matrix + backward_difference_matrix) / 2;
            end
            
            Forward_difference_cell{d}  = forward_difference_matrix;
            Backward_difference_cell{d} = backward_difference_matrix;
            Central_difference_cell{d}  = central_difference_matrix;
        end
        
    %% Matrix data selection function %%
        % This function works for matrices of arbitrary size
        function Matrix_Out = Matrix_Data_Selector(Matrix_In, dimension, number_dimensions, indices)
            dimension_selector = repmat({':'}, 1, number_dimensions);
            dimension_selector{dimension} = indices;

            Matrix_Out = Matrix_In(dimension_selector{:});
        end
end