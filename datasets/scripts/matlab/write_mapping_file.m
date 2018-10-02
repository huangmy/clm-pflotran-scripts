function write_mapping_file(filename,sid,did,wts,clm_nsoil_hydro,clm_nsoil_therm,clm_nsoil_mapped,pf_nsoil,pf_nsoil_mapped)

disp(filename)
fid=fopen(filename,'w');
fprintf(fid,'! Num of CLM soil layers (for soil hydrology)\n');
fprintf(fid,'clm_nlevsoi %d\n\n',clm_nsoil_hydro);
fprintf(fid,'! Num of CLM ground layers (for soil heat transport)\n');
fprintf(fid,'clm_nlevgrnd %d\n\n',clm_nsoil_therm);
fprintf(fid,'! Num of CLM layers mapped\n');
fprintf(fid,'clm_nlev_mapped %d\n\n',clm_nsoil_mapped);
fprintf(fid,'! Num of PFLOTRAN soil layers\n');
fprintf(fid,'pflotran_nlev %d\n\n',pf_nsoil);
fprintf(fid,'! Num of PFLOTRAN soil layers mapped\n');
fprintf(fid,'pflotran_nlev_mapped %d\n\n',pf_nsoil_mapped);
fprintf(fid,'! Num of weights\n');
fprintf(fid,'num_weights %d\n\n',size(sid,1));
fprintf(fid,'! FORMAT:\n');
fprintf(fid,'! pf_cell_idx clm_cell_idx weight\n');
fprintf(fid,'!\n');
fprintf(fid,'! Note: %spf_cell_idx%s and %sclm_cell_idx%s are in natural-order\n',char(39),char(39),char(39),char(39));
fprintf(fid,'!\n');for ii=1:size(sid,1)
    fprintf(fid,'%d %d %f\n',did(ii),sid(ii),wts(ii));
end
fclose(fid);
