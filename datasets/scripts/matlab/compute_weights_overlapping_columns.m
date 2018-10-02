function [source_cell_id, dest_cell_id, wts] = compute_weights_overlapping_columns(zv_dest, zv_source)



%
%    ----------          ----------
%
%                        ----------
%
%    ----------
%
%                        ----------
%
%    ----------
%
%                        ----------


count = 0;
source_cell_id = [];
dest_cell_id   = [];
wts            = [];

for ii = 1:length(zv_dest)-1
  
  loc = find(zv_source - zv_dest(ii)<=0);
  jj_beg = loc(end);
  
  loc = find(zv_source - zv_dest(ii+1)<=0);
  jj_end = loc(end)+1;
  
  jj_end = min(jj_end,length(zv_source));
   
  for jj = jj_beg:jj_end-1
    if (zv_source(jj)>= zv_dest(ii) && zv_source(jj+1)<=zv_dest(ii+1))
      weight = (zv_source(jj+1) - zv_source(jj))/(zv_dest(ii+1) - zv_dest(ii));
    elseif (zv_source(jj)< zv_dest(ii) && zv_source(jj+1)<=zv_dest(ii+1))
      weight = (zv_source(jj+1) - zv_dest(ii))/(zv_dest(ii+1) - zv_dest(ii));
    elseif (zv_source(jj)>= zv_dest(ii) && zv_source(jj+1)>zv_dest(ii+1))
      weight = (zv_dest(ii+1) - zv_source(jj))/(zv_dest(ii+1) - zv_dest(ii));
    elseif (zv_source(jj)< zv_dest(ii) && zv_source(jj+1)>zv_dest(ii+1))
      weight = 1.0;
    else
      error('Code failed!!!')
    end
    
    count = count + 1;
    dest_cell_id(count,1)    = ii;
    source_cell_id(count,1)  = jj;
    wts(count,1)             = weight;
  end
  
  
end
