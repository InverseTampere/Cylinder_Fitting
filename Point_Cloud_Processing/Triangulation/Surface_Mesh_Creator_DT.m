% This script creates a triangular surface mesh of a 3D point cloud using Delaunay triangulation

function [surface_triangle_mesh, num_facets] = Surface_Mesh_Creator_DT(point_cloud_matrix)

    %% Bounding mesh %%
        % 3D Delaunay triangulation
        Delaunay_Tetrahedrons   = delaunayTriangulation(point_cloud_matrix);     % This is a volumetric mesh, whilst only a mesh of the true surface is desired
        Tetrahedron_matrix      = Delaunay_Tetrahedrons.ConnectivityList;
        number_tetrahedrons     = size(Tetrahedron_matrix, 1);

        % The triangles within the tetrahedrons
        Triangle_fun    = @(tetra_ind) {nchoosek(tetra_ind, 3)};
        triangle_cell   = splitapply(Triangle_fun, Tetrahedron_matrix', 1:number_tetrahedrons);

        triangle_matrix_total   = vertcat(triangle_cell{:});

        % The number of times they occur (regardless of order of indices)
        triangle_matrix_sorted          = sort(triangle_matrix_total, 2);
        [triangle_matrix, ~, ind_list]  = unique(triangle_matrix_sorted, 'rows', 'stable');

        % Triangles on the surface are those that occur only once
        triangle_counts         = groupcounts(ind_list);
        surface_triangle_mesh   = triangle_matrix(triangle_counts == 1, :);        
        num_facets              = size(surface_triangle_mesh, 1);

end