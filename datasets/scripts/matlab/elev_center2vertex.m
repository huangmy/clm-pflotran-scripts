% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Given a 2D matrix of elevation at cell center, this function returns
% elevation values at vertices.
%
% Gautam Bisht (gbisht@lbl.gov)
% 01-02-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function zv_2d = elev_center2vertex(zz_2d)

[nx, ny] = size(zz_2d);

% allocate memory to save elevation values at cell vertex.
zv_2d(1:nx+1,1:ny+1) = 0;

% average elevation of 4 neighbors
zv_2d(2:nx,2:ny) = (zz_2d(1:end-1,1:end-1) + ...
    zz_2d(1:end-1,2:end  ) + ...
    zz_2d(2:end  ,1:end-1) + ...
    zz_2d(2:end  ,2:end  ))/4.0;

% average elevation of 2 neighbors
for ii = 1:nx-1
    jj = 1;
    zv_2d(ii+1,jj  ) = (zz_2d(ii,jj) + zz_2d(ii+1,jj))/2.0;
    jj = ny;
    zv_2d(ii+1,jj+1) = (zz_2d(ii,jj) + zz_2d(ii+1,jj))/2.0;
end

% average elevation of 2 neighbors
for jj = 1:ny-1
    ii = 1;
    zv_2d(ii  ,jj+1) = (zz_2d(ii,jj) + zz_2d(ii,jj+1))/2.0;
    ii = nx;
    zv_2d(ii+1,jj+1) = (zz_2d(ii,jj) + zz_2d(ii,jj+1))/2.0;
end

% set elevation at four corners to be equal to cell center
zv_2d(1   ,1   ) = zz_2d(1,1);
zv_2d(nx+1,1   ) = zz_2d(nx,1);
zv_2d(1   ,ny+1) = zz_2d(1,ny);
zv_2d(nx+1,ny+1) = zz_2d(nx,ny);


