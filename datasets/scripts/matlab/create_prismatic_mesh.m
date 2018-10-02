% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Returns connectivity information for an unstructured PFLOTRAN grid
% comprising of prismatic cells.
%
% Gautam Bisht (gbisht@lbl.gov)
% 07-08-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [vertices, cells, river, north, south, east, west] = ...
    create_prismatic_mesh( ...
    nx, ny, nz, ...
    xv_2d, yv_2d, zv_2d, ...
    dz_v, ...
    river_ii_beg, river_ii_end)

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Create prismatic PFLOTRAN mesh
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


xc_2d = (xv_2d(1:end-1,1:end-1)+xv_2d(2:end,1:end-1)+xv_2d(1:end-1,2:end)+xv_2d(2:end,2:end))/4;
yc_2d = (yv_2d(1:end-1,1:end-1)+yv_2d(2:end,1:end-1)+yv_2d(1:end-1,2:end)+yv_2d(2:end,2:end))/4;
zc_2d = (zv_2d(1:end-1,1:end-1)+zv_2d(2:end,1:end-1)+zv_2d(1:end-1,2:end)+zv_2d(2:end,2:end))/4;

% allocate memory to save (x,y,z) @ cell vertices
vertices(1:((nx+1)*(ny+1) + nx*ny)*(nz+1),1:3) = 0;

nv = 0;
for kk = 1:nz+1
    for jj = 1:ny+1
        for ii = 1:nx+1
            nv = nv+1;
            vertices(nv,1) = xv_2d(ii,jj);
            vertices(nv,2) = yv_2d(ii,jj);
            vertices(nv,3) = zv_2d(ii,jj) - dz_v(kk);
        end
    end
end

for kk = 1:nz+1
    for jj = 1:ny
        for ii = 1:nx
            nv = nv+1;
            vertices(nv,1) = xc_2d(ii,jj);
            vertices(nv,2) = yc_2d(ii,jj);
            vertices(nv,3) = zc_2d(ii,jj) - dz_v(kk);
        end
    end
end


ncells = 0;
nriver = 0;
neast  = 0;
nwest  = 0;
nnorth = 0;
nsouth = 0;

nv_cell = 6; % num of vertices forming a prismatic cell
nv_top  = 3; % num of vertices forming top face prismatic cell

nprism    = 4; % num of prisms within a hex

num_cells         = nx*ny*nz*nprism;         % Total number of prismatic cells within the domain
num_cells_per_lyr = num_cells/nz;            % Number of prismatic cells in each layer

% Cell IDs
layer_cell_ids(1:num_cells_per_lyr,1:nz)  = 0;
cells(1:num_cells,1:nv_cell+1)            = 0;
river(1:num_cells_per_lyr,1:nv_top+1)     = 0;
north(1:nx*nz,1:5)                        = 0;
south(1:nx*nz,1:5)                        = 0;
east( 1:ny*nz,1:5)                        = 0;
west( 1:ny*nz,1:5)                        = 0;


%  8--------7
%  |\     / |
%  | \   /  |
%  |  \ /   |
%  |   10   |
%  |   /\   |
%  |  /  \  |
%  | /    \ |
%  5--------6
%
%
%  4--------3
%  |\     / |
%  | \   /  |
%  |  \ /   |
%  |   9    |
%  |  / \   |
%  | /   \  |
%  |/     \ |
%  1--------2
%
%  Vertex order
%  1-2-9-5-6-10
%  2-3-9-6-7-10
%  3-4-9-7-8-10
%  4-1-9-8-5-10
%

for kk = 1:nz
    count_2d = 0;
    for jj = 1:ny
        for ii = 1:nx
            
            v1 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii;
            v2 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + 1;
            v3 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + 1 + (nx+1);
            v4 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + (nx+1);
            v5 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + (nx+1)*(ny+1);
            v6 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + 1 + (nx+1)*(ny+1);
            v7 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + 1 + (nx+1) + (nx+1)*(ny+1);
            v8 = (kk-1)*(nx+1)*(ny+1) + (jj-1)*(nx+1) + ii + (nx+1) + (nx+1)*(ny+1);
            
            v9 = (nx+1)*(ny+1)*(nz+1) + (kk-1)*(nx)*(ny) + (jj-1)*(nx) + ii;
            v10= (nx+1)*(ny+1)*(nz+1) + (kk-1)*(nx)*(ny) + (jj-1)*(nx) + ii + (nx)*(ny);
            
            
            ncells                       = ncells + 1;
            count_2d                     = count_2d + 1;
            layer_cell_ids(count_2d,kk)  = ncells;
            
            cells(ncells,1) = nv_cell;
            cells(ncells,2) = v1;
            cells(ncells,3) = v2;
            cells(ncells,4) = v9;
            cells(ncells,5) = v5;
            cells(ncells,6) = v6;
            cells(ncells,7) = v10;
            
            ncells                       = ncells + 1;
            count_2d                     = count_2d + 1;
            layer_cell_ids(count_2d,kk)  = ncells;
            
            cells(ncells,1) = nv_cell;
            cells(ncells,2) = v2;
            cells(ncells,3) = v3;
            cells(ncells,4) = v9;
            cells(ncells,5) = v6;
            cells(ncells,6) = v7;
            cells(ncells,7) = v10;
            
            
            ncells                       = ncells + 1;
            count_2d                     = count_2d + 1;
            layer_cell_ids(count_2d,kk)  = ncells;
            
            cells(ncells,1) = nv_cell;
            cells(ncells,2) = v3;
            cells(ncells,3) = v4;
            cells(ncells,4) = v9;
            cells(ncells,5) = v7;
            cells(ncells,6) = v8;
            cells(ncells,7) = v10;
            
            ncells                       = ncells + 1;
            count_2d                     = count_2d + 1;
            layer_cell_ids(count_2d,kk)  = ncells;
            
            cells(ncells,1) = nv_cell;
            cells(ncells,2) = v4;
            cells(ncells,3) = v1;
            cells(ncells,4) = v9;
            cells(ncells,5) = v8;
            cells(ncells,6) = v5;
            cells(ncells,7) = v10;
            
            if (ii == 1)
                nwest = nwest + 1;
                west(nwest,1)   = 4;
                west(nwest,2:5) = [v1 v5 v8 v9];
            end
            
            if (ii == nx)
                neast = neast + 1;
                east(neast,1)   = 4;
                east(neast,2:5) = [v2 v3 v7 v6];
            end
            
            if (jj == 1)
                nsouth = nsouth + 1;
                south(nsouth,1)   = 4;
                south(nsouth,2:5) = [v1 v2 v6 v5];
            end
            
            if (jj == ny)
                nnorth = nnorth + 1;
                north(nnorth,1)   = 4;
                north(nnorth,2:5) = [v3 v4 v8 v7];
            end
            
            if (kk == nz && ii >= river_ii_beg && ii <= river_ii_end)
                nriver = nriver + 1;
                river(nriver,:) = [3 v5 v6 v10];
                
                nriver = nriver + 1;
                river(nriver,:) = [3 v6 v7 v10];
                
                nriver = nriver + 1;
                river(nriver,:) = [3 v7 v8 v10];
                
                nriver = nriver + 1;
                river(nriver,:) = [3 v8 v5 v10];
            end
        end
    end
end
