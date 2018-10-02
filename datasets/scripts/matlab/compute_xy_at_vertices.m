% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Computes vertex coordinates (x,y) for a rectangular mesh with nx and ny
% grid points in x- and y-direction that are spaced at an uniform grid
% grid spacing of dx and dy.
%
% Gautam Bisht (gbisht@lbl.gov)
% 01-02-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [xv_2d,yv_2d] = compute_xy_at_vertices(dx,dy,nx,ny)

X = dx*nx;
Y = dx*ny;

xx = dx/2:dx:X;
yy = dy/2:dy:Y;

% allocate memory
xx_2d(1:nx  ,1:ny  ) = 0;
yy_2d(1:nx  ,1:ny  ) = 0;
xv_2d(1:nx+1,1:ny+1) = 0;
yv_2d(1:nx+1,1:ny+1) = 0;


for ii = 1:nx
    for jj = 1:ny
        % Coordinates @ cell centers
        xx_2d(ii,jj) = xx(ii);
        yy_2d(ii,jj) = yy(jj);
        
        % X-coordinates @ cell vertices
        xv_2d(ii  ,jj  ) = xx(ii)-dx/2;
        xv_2d(ii+1,jj  ) = xx(ii)+dx/2;
        xv_2d(ii  ,jj+1) = xx(ii)-dx/2;
        xv_2d(ii+1,jj+1) = xx(ii)+dx/2;
        
        % Y-coordinates @ cell vertices
        yv_2d(ii  ,jj  ) = yy(jj)-dy/2;
        yv_2d(ii+1,jj  ) = yy(jj)-dy/2;
        yv_2d(ii  ,jj+1) = yy(jj)+dy/2;
        yv_2d(ii+1,jj+1) = yy(jj)+dy/2;
    end
end

