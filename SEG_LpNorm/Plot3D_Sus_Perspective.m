% Create figure through 3D model
clear all
close all


addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB;

%% Input Files
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\Tile1';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Dipping_Plane';

meshfile = 'Mesh_20m.msh';

dsep = '\';

% obsfile = 'Tile_data.dat';
modelfile{1} = 'DippingPlane_Dip30_Azm30.sus';
modelfile{2} = 'MAG3D_TMI_lplq.sus';
modelfile{3} = 'MAG3D_TMI_l2l2.sus';

% mag_vecfile = '\l2l2\Mvec_TMVI_iter_.fld';
% predfile = '\l2l2\TMVI_iter_.pre';

% mag_vecfile = 'Tile1_MVI.fld';
% obsfile = 'Tile1_MAG3D_lplq.pre';
obsfile = 'MAG3D_TMI_l2l2.pre';

% Load surface
xyzS = load([work_dir dsep 'Plane_Surf.dat']);


% mag_vecfile = '..\..\m_rem.dat';
% obsfile = '..\Obs_loc_TMI.obs';

% mag_vecfile = '..\magvec.fld';
% obsfile = '..\..\Obs_RAW_REM_GRID_TMI.obs';
label2 ='$b^{TMI}\;(nT)$';

zpanel = 5;
ypanel = 8;

padE = 4;
padW = 4;

padN = 10;
padS = 12;

padT = 0;
padB = 10;

iso_cut_surf = 0.01;
iso_cut_vec = 0.01;

mmax = 0.2;

cam_ang = [45 15];

vscale = 1;

% Define cutting planes
% nvec(1,1:3) = [0 -0.2 1]; xo = 1000; yo = 850  ; zo = 1420;
% nvec(2,1:3) = [0 1 0]; xo(2) = 1000; yo(2) = 730  ; zo(2) = 1440;
nvec(3,1:3) = [1 0 0]; xo = 0; yo = 0  ; zo = 0;

nvec = spdiags( 1./ sqrt(sum(nvec.^2,2)) , 0 ,3 ,3) * nvec;

cut_name = ['A','B','C'];

% Color scheme

%% Load in model and plot
set(figure, 'Position', [200 50 750 400]); 

% Load data
[H, HI, HD, MI,MD, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir dsep obsfile]);
% H=1000;

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

% Move to local coordinates
obsx = obsx - median(xn);
obsy = obsy - median(yn);
obsz = obsz - zn(1);

xyzS(:,1) = xyzS(:,1) - median(xn);
xyzS(:,2) = xyzS(:,2) - median(yn);
xyzS(:,3) = xyzS(:,3) - zn(1);

tri = delaunay(xyzS(:,1),xyzS(:,2));


xn = xn - median(xn);
yn = yn - median(yn);
zn = zn - zn(1);

dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);
xx = (xn(2:end) + xn(1:end-1))/2;   xx = xx(padW+1:end-padE);
yy = (yn(2:end) + yn(1:end-1))/2;   yy = yy(padS+1:end-padN);
zz = (zn(2:end) + zn(1:end-1))/2;   zz = zz(padT+1:end-padB);

[XX,ZZ,YY] = meshgrid(xx,zz,yy); 

% Load magnetization vector
model = load([work_dir '\' modelfile{1}]);

m = reshape(model,nz,nx,ny);


mcell = size(m,1);
%% Load model 

m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

[YY2D,ZZ2D] = meshgrid(min(yy):20:max(yy),min(zz):20:max(zz)+min(dz));
XX2D = ones(size(YY2D)) * xo;

m2D = griddata(XX(:),YY(:),ZZ(:),m(:),XX2D,YY2D,ZZ2D,'nearest'); 
mmax = max(m2D(:));

ax1 = axes('Position',[0.0 .25 .6 .6]);

h = surf(XX2D,YY2D,ZZ2D,m2D); hold on
alpha(1)
view(cam_ang)
set(gca,'YDir','normal')

trisurf(tri,xyzS(:,1),xyzS(:,2),xyzS(:,3),'FaceColor','k','FaceAlpha',0.25);

% colormap('jet')
%plot_vec(ax1,XX2D,YY2D,ZZ2D,mx2D,my2D,mz2D,m2D,iso_cut_vec,0,0.5,1.25,3)

% axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])

set(gca,'YTickLabel',[],'Box','on')
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabh = get(gca,'ZLabel');

set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

set(gca, 'YAxisLocation', 'right')
set(gca,'YDir','normal')
hold on

[YY2D,ZZ2D] = meshgrid(yy+min(dx)/2,zz+min(dx)/2);
XX2D = ones(size(YY2D)) * xo;

cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear');
colormap(ax1,bb);
%colormap(ax2,bb);
caxis(ax1,[0 mmax])
% caxis(ax2,[0 mmax])
%% Add color bar
ax = axes('Position',[0.15 -0.1 .30 .30]);
cbar = colorbar(ax,'NorthOutside');

colormap(ax,bb);
caxis([0 mmax])

set(gca,'Visible','off');
text(0.5,2.2,'$\kappa\;(SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.45,1,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')


%% Plot Data
% Load data
ax2 = axes('Position',[0.475 .25 .6 .6]);

ndx = 100;
ndy = 100;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / ndx);
dy = (( ymax - ymin) / ndy);

x = xmin + cumsum(ones(1,ndx)*dx);
y = ymin + cumsum(ones(1,ndy)*dy);

[X,Y] = meshgrid(x,y);

data_interp     = griddata(obsx, obsy, d,X,Y,'linear'); 
% data_interp(isnan(data_interp)) = min(data_interp(:))*2;


h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(ax2,jet);
contour(x,y,data_interp,'k');
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
set(gca, 'XAxisLocation', 'top')
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[])
grid on
axis equal
axis tight
% title('$d^{Fwr}$', 'interpreter', 'latex','FontSize',14)
% text(min(xx)-dx*20, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)



%% Add colorbars
ax = axes('Position',[0.62 -0.1 .30 .30]);
cbar = colorbar('NorthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([min(data_interp(:)) max(data_interp(:))]))
set(gca,'Visible','off');
text(0.5,2.2,label2, 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.45,1,'$(b)$', 'interpreter', 'latex','FontSize',14)

%% Create second figure comparing recovered models
set(figure, 'Position', [200 50 750 400]); 

ax1 = axes('Position',[0.475 .25 .6 .6]);
% Load magnetization vector
model = load([work_dir '\' modelfile{2}]);

m = reshape(model,nz,nx,ny);

m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

[YY2D,ZZ2D] = meshgrid(min(yy):20:max(yy),min(zz):20:max(zz)+min(dz));
XX2D = ones(size(YY2D)) * xo;

m2D = griddata(XX(:),YY(:),ZZ(:),m(:),XX2D,YY2D,ZZ2D,'nearest'); 
mmax = max(m2D(:));

h = surf(XX2D,YY2D,ZZ2D,m2D); hold on
alpha(1)
view(cam_ang)
set(gca,'YDir','normal')

trisurf(tri,xyzS(:,1),xyzS(:,2),xyzS(:,3),'FaceColor','k','FaceAlpha',0.25);

% colormap('jet')
%plot_vec(ax1,XX2D,YY2D,ZZ2D,mx2D,my2D,mz2D,m2D,iso_cut_vec,0,0.5,1.25,3)

% axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])

set(gca,'YTickLabel',[],'Box','on')
% set(gca,'ZTickLabel',[],'Box','on')
% zlabel('$z$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabh = get(gca,'ZLabel');

% set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

set(gca, 'YAxisLocation', 'right')
set(gca,'YDir','normal')
hold on

[YY2D,ZZ2D] = meshgrid(yy+min(dx)/2,zz+min(dx)/2);
XX2D = ones(size(YY2D)) * xo;

cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear');
colormap(ax1,bb);
%colormap(ax2,bb);
caxis(ax1,[0 mmax])
% caxis(ax2,[0 mmax])

%%
ax2 = axes('Position',[0.0 .25 .6 .6]);

model = load([work_dir '\' modelfile{3}]);

m = reshape(model,nz,nx,ny);

m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );

[YY2D,ZZ2D] = meshgrid(min(yy):20:max(yy),min(zz):20:max(zz)+min(dz));
XX2D = ones(size(YY2D)) * xo;

m2D = griddata(XX(:),YY(:),ZZ(:),m(:),XX2D,YY2D,ZZ2D,'nearest'); 
mmax = max(m2D(:));

h = surf(XX2D,YY2D,ZZ2D,m2D); hold on
alpha(1)
view(cam_ang)
set(gca,'YDir','normal')

trisurf(tri,xyzS(:,1),xyzS(:,2),xyzS(:,3),'FaceColor','k','FaceAlpha',0.25);

% colormap('jet')
%plot_vec(ax1,XX2D,YY2D,ZZ2D,mx2D,my2D,mz2D,m2D,iso_cut_vec,0,0.5,1.25,3)

% axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])

set(gca,'YTickLabel',[],'Box','on')
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabh = get(gca,'ZLabel');

set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

set(gca, 'YAxisLocation', 'right')
set(gca,'YDir','normal')
hold on

[YY2D,ZZ2D] = meshgrid(yy+min(dx)/2,zz+min(dx)/2);
XX2D = ones(size(YY2D)) * xo;

cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear');
colormap(ax2,bb);
%colormap(ax2,bb);
caxis(ax2,[0 mmax])
% caxis(ax2,[0 mmax])
%% Add color bar
ax = axes('Position',[0.35 -0.1 .30 .30]);
cbar = colorbar(ax,'NorthOutside');

colormap(ax,bb);
caxis([0 mmax])

set(gca,'Visible','off');
text(0.5,2.2,'$\kappa\;(SI)$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(-.25,1.75,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
text(1.5,1.75,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')