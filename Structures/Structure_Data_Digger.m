% This script goes through the structure to retrieve and aggregate the data into an array
% The data is given in a single-level output structure as well as an array for which each variable spans the indicated columns

function [Data_Structure, data_fields, number_datasets, data_array, Index_Structure] = Structure_Data_Digger(Input_Structure)

    %% Unfold the input structure %%
        % This results in a single-level structure
        Data_Structure  = Data_Extraction(Input_Structure);
        data_fields     = fieldnames(Data_Structure);
        number_datasets = length(data_fields);

    %% Convert to array format %%
        % The structure is converted into an array
        data_cell   = struct2cell(Data_Structure);
        data_array  = horzcat(data_cell{:});

        % The indices corresponding to each variable are saved with their names in a structure
        Dim_fun         = @(array) size(array, 2);          % Number of dimensions of a variable equals its number of columns
        num_dim_list    = cellfun(Dim_fun, data_cell);
        cumsum_dim_list = [0; cumsum(num_dim_list)];

        Index_Structure = struct();

        for d = 1 : number_datasets
            % The name of this data
            data_field = data_fields{d};

            % The columns it inhibits
            first_column    = cumsum_dim_list(d) + 1;
            last_column     = cumsum_dim_list(d + 1);
            columns         = first_column : last_column;

            % Saved in the structure
            Index_Structure.(data_field) = columns;
        end

    %% Extraction function %%
        % Generates a structure containing the data of each variable 
    function Output_Structure = Data_Extraction(Input_Structure)
            % Number of elements of the structure
            number_elements = numel(Input_Structure);
            
            % Loop through the elements to gather the data
            Output_Structure = struct();
        
            for i = 1 : number_elements
                % Fields of the i'th structure
                field_names     = fieldnames(Input_Structure(i));
                number_fields   = length(field_names);
            
                for f = 1 : number_fields
                    % This field's data
                    field_name = field_names{f};
                    field_data = Input_Structure(i).(field_name);
        
                    % If it is a structure, the unfolding script is called recursively at a deeper level
                    if isstruct(field_data)
                        % Output structure of the recursively called function
                        Output_Structure_r  = Data_Extraction(field_data);
                        Output_Structure    = Structure_Merging(Output_Structure, Output_Structure_r);
        
                    % Otherwise if it is a cell or numeric array, an output is generated
                    elseif isnumeric(field_data) || iscell(field_data)
                        % If the field already exists, the data is added to it
                        if isfield(Output_Structure, field_name) 
                            existing_field_data = Output_Structure.(field_name);
                        else
                            existing_field_data = [];
                        end
                        
                        Output_Structure.(field_name) = [existing_field_data; field_data];
        
                    else
                        class_type = class(field_data);
                        error('The structure''s data is of type %s. It must be either cell or numeric.', class_type);
                    end
                end
            end
        end
    
    %% Merging function %%
        % Merges the data of two single-level structures
        function Merged_Structure = Structure_Merging(Structure_1, Structure_2)
        % Fields in both structures
        field_names_2   = fieldnames(Structure_2);
        number_fields_2 = length(field_names_2);
    
        % Cycling through the fields
        Merged_Structure = Structure_1;
    
        for f = 1 : number_fields_2
            field_name = field_names_2{f};
            field_data = Structure_2.(field_name);
    
            % If the field already exists in the first structure, the data is appended
            if isfield(Structure_1, field_name)
                existing_data   = Merged_Structure.(field_name);
                merged_data     = [existing_data; field_data];
            
                Merged_Structure.(field_name) = merged_data;
    
            % Otherwise a new field is added
            else
                Merged_Structure.(field_name) = field_data;
            end
        end
    end
end