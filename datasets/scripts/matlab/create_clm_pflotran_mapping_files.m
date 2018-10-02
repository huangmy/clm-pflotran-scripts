function create_clm_pflotran_mapping_files(resolution_txt)

% Determin the resolution of the problem
switch resolution_txt
  case '20m'
    nx_clm = 20;
    ny_clm = 20;
    nz_clm = 10;
    
    nx_pf = 20;
    ny_pf = 20;
    nz_pf = 62;
    
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
    
    dx_pf = 2.0; % [m]
    dy_pf = 2.0; % [m]
    dz_pf = 0.5; % [m]
    
    pflotran_material_hdf5_fname = '../../pflotran_2m/clmpf_400x400x31_2mx2mxhalfRes_material_mapped.h5';
  otherwise
    error(['Unsupported resolution: ' resolution_txt])
end



% Do some error checking
if ((nx_clm ~= nx_pf))
  error('The code cannot account for unequal number of grid cells in X-direction between CLM and PFLOTRAN')
end

if ((ny_clm ~= ny_pf))
  error('The code cannot account for unequal number of grid cells in Y-direction between CLM and PFLOTRAN')
end


% Determine dz for CLM grid
dzs_clm = get_clm_dzs (10)';

% Determine dz for PFLOTRAN grid
dzs_pf  = ones(nz_pf,1)*dz_pf;

% The z of cell edges
zvert_clm = cumsum([0; dzs_clm]);
zvert_pf  = cumsum([0; dzs_pf ]);

[col_p2c_sid, col_p2c_did, col_p2c_wts] = compute_weights_overlapping_columns(zvert_clm, zvert_pf );
[col_c2p_sid, col_c2p_did, col_c2p_wts] = compute_weights_overlapping_columns(zvert_pf , zvert_clm);



pf_cell_ids = zeros(nx_pf,ny_pf,nz_pf);
clm_cell_ids = zeros(nx_clm,ny_clm,nz_clm);

count = 0;
for kk = 1:nz_pf
  for jj = 1:ny_pf
    for ii = 1:nx_pf
      count = count + 1;
      pf_cell_ids(ii,jj,kk) = count;
    end
  end
end


count = 0;
for jj = 1:ny_clm
  for ii = 1:nx_clm
    for kk = nz_clm:-1:1
      count = count + 1;
      clm_cell_ids(ii,jj,kk) = count;
    end
  end
end

% Read in PFLOTRAN material dataset
pflotran_dataset_cell_id = double(h5read(pflotran_material_hdf5_fname,'/Materials/Cell Ids'));
pflotran_dataset_mat_id = double(h5read(pflotran_material_hdf5_fname,'/Materials/Material Ids'));

% perform a check
if (nx_pf*ny_pf*nz_pf ~= length(pflotran_dataset_cell_id))
  error('The length of PFLOTRAN material dataset does not match nxXnyXnz');
end

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


%
p2c_sid = zeros(length(col_p2c_did)*nx_clm*ny_clm,1);
p2c_did = zeros(length(col_p2c_did)*nx_clm*ny_clm,1);
p2c_wts = zeros(length(col_p2c_did)*nx_clm*ny_clm,1);

c2p_sid = zeros(length(col_c2p_did)*nx_clm*ny_clm,1);
c2p_did = zeros(length(col_c2p_did)*nx_clm*ny_clm,1);
c2p_wts = zeros(length(col_c2p_did)*nx_clm*ny_clm,1);


% SUBSURFACE: PFLOTRAN -----> CLM
nsid = length(col_p2c_sid);
ndid = length(col_p2c_did);

count = 0;
for jj = 1:ny_pf
  for ii = 1:nx_pf
    
    kk = pf_surface_cell_kindex(ii,jj);
    sid = reshape(pf_cell_ids(ii,jj,kk - col_p2c_sid + 1),nsid,1);
    did = reshape(clm_cell_ids(ii,jj,nz_clm - col_p2c_did + 1),ndid,1);
    
    count = count + 1;
    p2c_sid((count-1)*nsid+1:count*nsid) = sid;
    p2c_did((count-1)*nsid+1:count*nsid) = did;
    p2c_wts((count-1)*nsid+1:count*nsid) = col_p2c_wts;
  end
end

% SUBSURFACE:  CLM -----> PFLOTRAN
nsid = length(col_c2p_sid);
ndid = length(col_c2p_did);

count = 0;
for jj = 1:ny_pf
  for ii = 1:nx_pf
    
    kk = pf_surface_cell_kindex(ii,jj);
    sid = reshape(clm_cell_ids(ii,jj,nz_clm - col_c2p_did + 1),ndid,1);
    did = reshape(pf_cell_ids(ii,jj,kk - col_c2p_sid + 1),nsid,1);
    
    count = count + 1;
    c2p_sid((count-1)*nsid+1:count*nsid) = sid;
    c2p_did((count-1)*nsid+1:count*nsid) = did;
    c2p_wts((count-1)*nsid+1:count*nsid) = col_c2p_wts;
  end
end

system('mkdir -p mapping_files');

filename = sprintf('mapping_files/pf2clm_map_%dx%dx%dPFLOTRAN_%dx%dx%dCLM_%s.meshmap',nx_pf,ny_pf,nz_pf,nx_clm,ny_clm,nz_clm,datestr(now, 'cyymmdd'));
[~,idx] = sort(p2c_did);
p2c_sid = p2c_sid(idx);
p2c_did = p2c_did(idx);
p2c_wts = p2c_wts(idx);
write_mapping_file(filename,p2c_sid,p2c_did,p2c_wts,nz_clm,nz_clm,nz_clm,nz_pf,nz_pf)

filename = sprintf('mapping_files/clm2pf_map_%dx%dx%dPFLOTRAN_%dx%dx%dCLM_%s.meshmap',nx_pf,ny_pf,nz_pf,nx_clm,ny_clm,nz_clm,datestr(now, 'cyymmdd'));
[~,idx] = sort(c2p_did);
c2p_sid = c2p_sid(idx);
c2p_did = c2p_did(idx);
c2p_wts = c2p_wts(idx);
write_mapping_file(filename,c2p_sid,c2p_did,c2p_wts,nz_clm,nz_clm,nz_clm,nz_pf,nz_pf)


% SUBSURFACE:  CLM -----> PFLOTRAN
surf_c2p_sid = reshape(clm_cell_ids(:,:,end),nx_clm*ny_clm,1);
surf_c2p_did = reshape(pf_surface_cell_ids,nx_clm*ny_clm,1);
surf_c2p_wts = ones(nx_clm*ny_clm,1);

filename = sprintf('mapping_files/clm2pf_map_surface_%dx%dx%dPFLOTRAN_%dx%dx%dCLM_%s.meshmap',nx_pf,ny_pf,nz_pf,nx_clm,ny_clm,nz_clm,datestr(now, 'cyymmdd'));
write_mapping_file(filename,surf_c2p_sid,surf_c2p_did,surf_c2p_wts,nz_clm,nz_clm,nz_clm,nz_pf,nz_pf)


