% The given cell array of booleans is applied to the rows of the given cell array containing matrices

function filtered_data_cell = Cell_Array_Boolean(data_cell, bool_cell)

    % Filtering
    Filtering_fun       = @(data_matrix, bool_list) data_matrix(bool_list, :);
    filtered_data_cell  = cellfun(Filtering_fun, data_cell, bool_cell, 'UniformOutput', false);

end