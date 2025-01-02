% Edge vertices are determined of the given mesh which need not be triangular

function [edge_vertex_matrix, number_edge_vertices] = Mesh_Edge_Vertices(mesh)
    
    %% Vertices %%
        % The vertices within the mesh        
        Vertex_Ind_fun  = @(polygon_ind) {nchoosek(polygon_ind, 2)};
        number_polygons = size(mesh, 1);
        vertex_ind_cell = splitapply(Vertex_Ind_fun, mesh', 1 : number_polygons);

    %% Edge vertices %%
        % The number of times each vertex occurs (regardless of order of indices)
        vertex_ind_matrix                   = vertcat(vertex_ind_cell{:});
        vertex_ind_matrix                   = sort(vertex_ind_matrix, 2);
        [vertex_ind_matrix, ~, ind_list]    = unique(vertex_ind_matrix, 'rows', 'stable');

        % Edge vertices on the surface are those that occur only once
        vertex_counts           = groupcounts(ind_list);
        edge_vertex_matrix      = vertex_ind_matrix(vertex_counts == 1, :);        
        number_edge_vertices    = size(edge_vertex_matrix, 1);

end