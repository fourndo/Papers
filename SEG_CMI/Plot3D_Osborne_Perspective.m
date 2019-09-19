% Create figure through 3D model
clear all



addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB;

%% Input Files
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\Tile1';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock\CMI\MVI';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Osborne\Inversion\ROT40\CMI\MVI';

meshfile = 'Tile1.msh';

dsep = '\';

% obsfile = 'Tile_data.dat';
% model_true = '..\..\Block.sus';
geo_unit = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Osborne\Section_21360.dat';

% mag_vecfile = '\l2l2\Mvec_TMVI_iter_.fld';
% predfile = '\l2l2\TMVI_iter_.pre';
mod_file = 'Tile1_MAG3D_l2l2.sus';
mag_vecfile = 'Tile1_MVI.fld';
obsfile = 'Tile1_MVI.pre';

% mag_vecfile = '..\..\m_rem.dat';
% obsfile = '..\Obs_loc_TMI.obs';

% mag_vecfile = '..\magvec.fld';
% obsfile = '..\..\Obs_RAW_REM_GRID_TMI.obs';


zpanel = 5;
ypanel = 8;

padE = 4;
padW = 13;

padN = 12;
padS = 12;

padT = 0;
padB = 6;

x_ext = [600  1350];

iso_cut_surf = 0.0075;
iso_cut_vec = 0.0075;

mmax = 0.5;

cam_ang = [0 0];

vscale = 1;

% Define cutting planes
% nvec(1,1:3) = [0 -0.2 1]; xo = 1000; yo = 850  ; zo = 1420;
% nvec(2,1:3) = [0 1 0]; xo(2) = 1000; yo(2) = 730  ; zo(2) = 1440;
yplane = 30;

nvec(3,1:3) = [0 1 0]; 

nvec = spdiags( 1./ sqrt(sum(nvec.^2,2)) , 0 ,3 ,3) * nvec;

cut_name = ['A','B','C'];

% Color scheme

%% Load in model and plot
set(figure, 'Position', [200 50 750 400]); 

ax1 = axes('Position',[0.3 .05 .9 .9]);
% Load geo unit
bif = load(geo_unit);

% Load data
[H, HI, HD, MI,MD, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir dsep obsfile]);
% H=1000;

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

% Move to local coordinates
bif(:,1) = bif(:,1) - xn(1);
bif(:,2) = bif(:,2) - yn(1);
bif(:,3) = bif(:,3) - zn(1);

obsx = obsx - xn(1);
obsy = obsy - yn(1);
obsz = obsz - zn(1);

x0 = 455970; y0 = yn(yplane)  ; z0 = 300;


x0 = x0 - xn(1);
y0 = y0 - yn(1);
z0 = z0 - zn(1);

xn = xn - xn(1); 
yn = yn - yn(1); 
zn = zn - zn(1); 

xx = (xn(2:end) + xn(1:end-1))/2;   xx = xx(padW+1:end-padE);
yy = (yn(2:end) + yn(1:end-1))/2;   yy = yy(padS+1:end-padN);
zz = (zn(2:end) + zn(1:end-1))/2;   zz = zz(padT+1:end-padB);

[XX,ZZ,YY] = meshgrid(xx,zz,yy); 

% Load mag true
% m_true = load([work_dir '\' model_true]);
% m_true = reshape(m_true,nz,nx,ny);


% Load magnetization vector
mag_model = load([work_dir '\' mag_vecfile]);


% Load inverted model
% m = load([work_dir '\..\..\Effec_sus_20mGrid.sus']);
m = sqrt(sum(mag_model.^2,2));
m = reshape(m,nz,nx,ny);


mcell = size(mag_model,1);
%% Create model magnetization vectors
% Azimuth and dip of magnitization
if size(mag_model,2)==2
    
    mag_xyz = azmdip_2_xyz( mag_model(:,1)+180 , mag_model(:,2) );
    

else
    
    mag_xyz = mag_model;
%     mag_xyz(mag_xyz(:,3)~=0,1:2)=0;
%     mag_xyz(mag_xyz(:,3)~=0,3)=m(mag_xyz(:,3)~=0);
%     M = [spdiags(m_true(:).*mag_model(:,1),0,mcell,mcell);spdiags(m_true(:).*mag_model(:,2),0,mcell,mcell);spdiags(m_true(:).*mag_model(:,3),0,mcell,mcell)];

   
end

mx = reshape(mag_xyz(:,1),nz,nx,ny);
my = reshape(mag_xyz(:,2),nz,nx,ny);
mz = reshape(-mag_xyz(:,3),nz,nx,ny);

mx = mx( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
my = my( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
mz = mz( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
    
% m_true = m_true( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );




%%



[XX2D,ZZ2D] = meshgrid(min(xx):10:max(xx),min(zz):10:max(zz));
YY2D = ones(size(XX2D)) * y0;
% Interpolate on cutting plane
% m2D = get_model_top(m1,nx,ny,nz,-100);

% cvec = [min(m(:)) prctile(m(m>0),[10 30 60]) max(m(:))];
% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0]/255,sort(m(:)));

m2D = griddata(XX(:),YY(:),ZZ(:),m(:),XX2D,YY2D,ZZ2D,'natural'); 
% m2D = F(XX2D,YY2D,ZZ2D);

% F = TriScatteredInterp(XX(:),YY(:),ZZ(:),mx(:),'natural'); 
mx2D = squeeze(mx(:,:,yplane-padS));%F(XX2D,YY2D,ZZ2D);

% F = TriScatteredInterp(XX(:),YY(:),ZZ(:),my(:),'natural'); 
my2D = squeeze(my(:,:,yplane-padS));%F(XX2D,YY2D,ZZ2D);

% F = TriScatteredInterp(XX(:),YY(:),ZZ(:),mz(:),'natural'); 
mz2D = squeeze(mz(:,:,yplane-padS));%F(XX2D,YY2D,ZZ2D);

% F = TriScatteredInterp(XX(:),YY(:),ZZ(:),m_true(:),'natural'); 
% m2D_true = F(XX2D,YY2D,ZZ2D);



% temp = (temp');
vec = (mx2D.^2+my2D.^2+mz2D.^2).^0.5;

% Scale large vectors to mmax
indx = vec > mmax;
scal = 1 / max(vec(indx));

mx2D(indx) = mx2D(indx)*scal; 
my2D(indx) = my2D(indx)*scal;
mz2D(indx) = mz2D(indx)*scal;
vec(indx) = vec(indx) *scal;
vec = vec(:);

h = surf(XX2D,YY2D,ZZ2D,m2D,'EdgeColor','none'); hold on
alpha(0.85)
view(cam_ang)
set(gca,'YDir','normal')



% set(gca,'YTickLabel',[],'Box','on')
% set(gca,'XTickLabel',[],'Box','on')
% zlabel('$z$', 'interpreter', 'latex','FontSize',14)
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
set(gca,'ZTickLabel',[])
% zlabh = get(gca,'ZLabel');
% set(zlabh,'Position',get(zlabh,'Position') + [50 0 0])

% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [1 0 0])

% set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

% ax1.XTick = [x_ext(1)+25 x_ext(2)-25] ;
% ax1.XTickLabel = {num2str(x_ext(1)+25) ,num2str(x_ext(2)-25)};
set(gca, 'XAxisLocation', 'top')
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(gca,'YDir','normal')
hold on

[XX2D,ZZ2D] = meshgrid(xx+min(dx)/2,zz+min(dx)/2);
YY2D = ones(size(XX2D)) * y0;
arrow3([XX2D(:),YY2D(:),ZZ2D(:)],[mx2D(:)./vec,my2D(:)./vec,mz2D(:)./vec],'k',vec*0.005/max(vec),vec*0.01/max(vec),'cone')

s = scatter3(bif(:,1),bif(:,2)-5,bif(:,3),'k.','MarkerEdgeColor','k');

% Equalize color
[nn,hx] = histcounts(m2D(:),5);
y = round( ( nn - min(nn) ) / (length(m2D(:)) - min(nn) ) * 4 ) + 1;
cvec = mmax*[0 cumsum(y)/sum(y)*1.05];
bb = interp1(cvec',[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear','extrap');
colormap(ax1,bb);
%colormap(ax2,bb);
caxis(ax1,[0 mmax])
ylim([y0-10, y0 + 10])
zlim([min(zz), max(zz)])
xlim(x_ext)
% caxis(ax2,[0 mmax])
%% Add color bar

ax = axes('Position',[0.375 -0.05 .30 .30]);
cbar = colorbar(ax,'NorthOutside');

colormap(ax,bb);
caxis([0 mmax])

set(gca,'Visible','off');
text(0.5,2.2,'$\kappa_{e}$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
text(-.4,1.75,'$(a)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')
text(1.5,1.75,'$(b)$', 'interpreter', 'latex','FontSize',14,'HorizontalAlignment','center')

%% Plot k effective
ax2 = axes('Position',[0.05 .25 .5 .5]);

model = load([work_dir '\' mod_file]);

m = reshape(model,nz,nx,ny);

m = m( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
%m(m==-100)=nan;

[XX2D,ZZ2D] = meshgrid(min(xx):10:max(xx),min(zz):10:max(zz));
YY2D = ones(size(XX2D)) * y0;

m2D = griddata(XX(:),YY(:),ZZ(:),m(:),XX2D,YY2D,ZZ2D,'natural'); 
% mmax = 0.5;

h = surf(XX2D,YY2D,ZZ2D,m2D,'EdgeColor','none'); hold on
alpha(1)
view(cam_ang)
set(gca,'YDir','normal')
caxis(ax2,[0 mmax])
% trisurf(tri,xyzS(:,1),xyzS(:,2),xyzS(:,3),'FaceColor','k','FaceAlpha',0.25);

% colormap('jet')
%plot_vec(ax1,XX2D,YY2D,ZZ2D,mx2D,my2D,mz2D,m2D,iso_cut_vec,0,0.5,1.25,3)
% cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
% bb = interp1(cvec',[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-2:mmax,'linear','extrap');
% axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])

set(gca,'YTickLabel',[],'Box','on')
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabh = get(gca,'ZLabel');

set(get(gca,'ZLabel'),'Rotation',360);
axis equal
% grid on
% 
set(gca, 'XAxisLocation', 'top')
% set(gca,'YDir','normal')
% hold on

% [YY2D,ZZ2D] = meshgrid(yy+min(dx)/2,zz+min(dx)/2);
% XX2D = ones(size(YY2D)) * xo;
scatter3(bif(:,1),bif(:,2),bif(:,3),'k.','MarkerEdgeColor','k')



colormap(ax2,bb);
ylim([y0-10, y0 + 10])
zlim([min(zz), max(zz)])
xlim(x_ext)
% caxis(ax2,[0 mmax])
