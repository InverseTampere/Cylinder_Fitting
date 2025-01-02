% This script merges two structures with different fields containing numerical data
% Note that the structure vertically concatenates the given data into matrix format

function [merged_structure, merged_fieldnames, number_fields] = Different_Structure_Merger(structure_a, structure_b)
    
    %% Merging the two structures %%
        % The fieldnames of both structures
        merged_fieldnames   = [fieldnames(structure_a); fieldnames(structure_b)];
        number_fields       = length(merged_fieldnames);
        
        % Cell array of the contained data
        Merging_Function_handle = @(fieldname) Merging_Function(fieldname);
        structure_data_cell     = cellfun(Merging_Function_handle, merged_fieldnames, 'UniformOutput', false);
        
        % The new structure
        merged_structure    = cell2struct(structure_data_cell, merged_fieldnames, 1);
        
    %% Merging function %%
    function data_matrix = Merging_Function(fieldname)
        % It is attempted to retrieve the data from the first structure
        try
            data_cell = {structure_a.(fieldname)};
        
        % If that does not succeed, it means the data is contained in the second structure
        catch
            data_cell = {structure_b.(fieldname)};
        end
        
        % The cell array is converted to matrix format
        data_matrix = vertcat(data_cell{:});
    end
end