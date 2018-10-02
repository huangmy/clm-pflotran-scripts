function create_ugrid(resolution_txt, mesh_compatibility_format)

CLM_PFLTORAN_COMPATIABLE_MESH = 1;
PFLTORAN_DEV_COMPATIABLE_MESH = 2;

% Determin the resolution of the problem
switch resolution_txt
  case '20m'
    nx_clm = 20;
    ny_clm = 20;
    nz_clm = 10;
    
    nx_pf = 20;
    ny_pf = 20;
    nz_pf = 62;

    river_ii_beg = 15;
    river_ii_end = 20;
    
    origin_x = 0;
    origin_y = 0;
    origin_z = 90;
    
    dx_pf = 20.0; % [m]
    dy_pf = 20.0; % [m]
    dz_pf = 0.5; % [m]
    
    pflotran_material_hdf5_fname = '../../pflotran_20m/clmpf_400x400x31_20mx20mxhalfRes_material_mapped.h5';
  case '10m'
    nx_clm = 40;
    ny_clm = 40;
    nz_clm = 10;
    
    nx_pf = 40;
    ny_pf = 40;
    nz_pf = 62;
    
    river_ii_beg = 30;
    river_ii_end = 40;
    
    origin_x = 0;
    origin_y = 0;
    origin_z = 90;

    dx_pf = 10.0; % [m]
    dy_pf = 10.0; % [m]
    dz_pf = 0.5; % [m]
    
    pflotran_material_hdf5_fname = '../../pflotran_10m/clmpf_400x400x31_10mx10mxhalfRes_material_mapped.h5';
  case '2m'
    nx_clm = 200;
    ny_clm = 200;
    nz_clm = 10;
    
    nx_pf = 200;
    ny_pf = 200;
    nz_pf = 62;
    
    river_ii_beg = 150;
    river_ii_end = 200;
    
    origin_x = 0;
    origin_y = 0;
    origin_z = 90;

    dx_pf = 2.0; % [m]
    dy_pf = 2.0; % [m]
    dz_pf = 0.5; % [m]
    
    pflotran_material_hdf5_fname = '../../pflotran_2m/clmpf_400x400x31_2mx2mxhalfRes_material_mapped.h5';
  otherwise
    error(['Unsupported resolution: ' resolution_txt])
end

% Determine the HDF5 is for which model: clm-pflotran or pflotran-dev?
% This is necessary because of the recent changes in pflotran-dev.
switch lower(mesh_compatibility_format)
    case 'clm_pflotran'
        mesh_format = CLM_PFLTORAN_COMPATIABLE_MESH;
    case 'pflotran_dev'
        mesh_format = PFLTORAN_DEV_COMPATIABLE_MESH;
    otherwise
        error('Invalid entries for mesh_compatibility_format.');
end


% Read in PFLOTRAN material dataset
pflotran_dataset_cell_id = double(h5read(pflotran_material_hdf5_fname,'/Materials/Cell Ids'));
pflotran_dataset_mat_id = double(h5read(pflotran_material_hdf5_fname,'/Materials/Material Ids'));

% perform a check
if (nx_pf*ny_pf*nz_pf ~= length(pflotran_dataset_cell_id))
  error('The length of PFLOTRAN material dataset does not match nxXnyXnz');
end

pf_cell_ids = reshape([1:nx_pf*ny_pf*nz_pf],nx_pf,ny_pf,nz_pf);

% Determine which PFLOTRAN cells are active;
pf_active_cells = ones(nx_pf,ny_pf,nz_pf);
loc = find (pflotran_dataset_mat_id<= 0);
pf_active_cells(pflotran_dataset_cell_id(loc)) = 0;


% Find IDs of active cells on the surface of PFLOTRAN grid.
pf_surface_cell_ids    = zeros(nx_pf,ny_pf);
pf_surface_cell_kindex = zeros(nx_pf,ny_pf);
for jj = 1:ny_pf
  for ii = 1:nx_pf
    loc = find(pf_active_cells(ii,jj,:) == 1);
    pf_surface_cell_ids(ii,jj)    = pf_cell_ids(ii,jj,loc(end));
    pf_surface_cell_kindex(ii,jj) = loc(end);
  end
end


zc_2d = pf_surface_cell_kindex*dz_pf;
zv_2d = elev_center2vertex(zc_2d);

[xv_2d,yv_2d] = compute_xy_at_vertices(dx_pf,dy_pf,nx_pf,ny_pf);

dz_v = flipud(cumsum(ones(nz_pf+1,1)*0.5));

[vertices, cells, river, north, south, east, west] = create_prismatic_mesh(nx_pf,ny_pf,nz_pf,xv_2d,yv_2d,zv_2d,dz_v,river_ii_beg,river_ii_end);

switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        suffix = 'clm_pf';
    case PFLTORAN_DEV_COMPATIABLE_MESH
        suffix = 'pf_dev';
end

filename = sprintf('ugrid/pflotran_ugrid_%04.2fmx%04.2fmx%04.2fm_%s_%s',dx_pf,dy_pf,dz_pf,suffix,datestr(now, 'cyymmdd'));
loc = find(filename == '.');
filename(loc) = 'p';
filename = [filename '.h5'];

disp(filename)

system(['rm -f ' filename]);

vertices(:,1) = vertices(:,1) + origin_x;
vertices(:,2) = vertices(:,2) + origin_y;
vertices(:,3) = vertices(:,3) + origin_z;

% Regions
ncells = size(cells,1);
nv = size(vertices,1);

all      = [1:ncells];
% top_edge = [south;east;north;west];

% Save data in HDF5 file

% cells
h5create(filename,'/Domain/Cells',[size(cells,2) ncells],'Datatype','uint64');
h5write(filename,'/Domain/Cells',uint64(cells)',[1 1],[size(cells,2) ncells]);

% vertices
h5create(filename,'/Domain/Vertices',[3 nv]);
h5write(filename,'/Domain/Vertices',(vertices'),[1,1],[3 nv]);

% region: all
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/all',[ncells],'Datatype','uint64');
        h5write(filename,'/Regions/all',uint64(all),1,[ncells]);
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/all/Cell Ids',[ncells],'Datatype','uint64');
        h5write(filename,'/Regions/all/Cell Ids',uint64(all),1,[ncells]);
end

% % region: top
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/river_on_top_surface',fliplr(size(river)),'Datatype','uint64');
        h5write(filename,'/Regions/river_on_top_surface',uint64(river'),[1,1],fliplr(size(river)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/river_on_top_surface/Vertex Ids',fliplr(size(river)),'Datatype','uint64');
        h5write(filename,'/Regions/river_on_top_surface/Vertex Ids',uint64(river'),[1,1],fliplr(size(river)));
end


% % region: north
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/north',fliplr(size(north)),'Datatype','uint64');
        h5write(filename,'/Regions/north',uint64(north'),[1,1],fliplr(size(north)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/north/Vertex Ids',fliplr(size(north)),'Datatype','uint64');
        h5write(filename,'/Regions/north/Vertex Ids',uint64(north'),[1,1],fliplr(size(north)));
end


% % region: south
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/south',fliplr(size(south)),'Datatype','uint64');
        h5write(filename,'/Regions/south',uint64(south'),[1,1],fliplr(size(south)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/south/Vertex Ids',fliplr(size(south)),'Datatype','uint64');
        h5write(filename,'/Regions/south/Vertex Ids',uint64(south'),[1,1],fliplr(size(south)));
end

% % region: west
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/west',fliplr(size(west)),'Datatype','uint64');
        h5write(filename,'/Regions/west',uint64(west'),[1,1],fliplr(size(west)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/west/Vertex Ids',fliplr(size(west)),'Datatype','uint64');
        h5write(filename,'/Regions/west/Vertex Ids',uint64(west'),[1,1],fliplr(size(west)));
end

% % region: east
switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/east',fliplr(size(east)),'Datatype','uint64');
        h5write(filename,'/Regions/east',uint64(east'),[1,1],fliplr(size(east)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/east/Vertex Ids',fliplr(size(east)),'Datatype','uint64');
        h5write(filename,'/Regions/east/Vertex Ids',uint64(east'),[1,1],fliplr(size(east)));
end


switch mesh_format
    case CLM_PFLTORAN_COMPATIABLE_MESH
        h5create(filename,'/Regions/river_on_east',fliplr(size(east)),'Datatype','uint64');
        h5write(filename,'/Regions/river_on_east',uint64(east'),[1,1],fliplr(size(east)));
    case PFLTORAN_DEV_COMPATIABLE_MESH
        h5create(filename,'/Regions/river_on_east/Vertex Ids',fliplr(size(east)),'Datatype','uint64');
        h5write(filename,'/Regions/river_on_east/Vertex Ids',uint64(east'),[1,1],fliplr(size(east)));
end



% write global attributes
[~,user_name]=system('echo $USER');
h5writeatt(filename,'/','Created_by: ' ,user_name(1:end-1));
h5writeatt(filename,'/','Created_on: ' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '));



