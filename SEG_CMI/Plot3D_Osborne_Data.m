% Create figure through 3D model
clear all



addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB';

%% Input Files
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_AMI\Tile1';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1';
% work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\SingleBlock\CMI\MVI';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Osborne\Inversion\UTM';

meshfile = 'Mesh_20m_ROT40.msh';

dsep = '\';

% obsfile = 'Tile_data.dat';
% model_true = '..\..\Block.sus';
geo_unit = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Osborne\Section_21360_UTM.dat';

% mag_vecfile = '\l2l2\Mvec_TMVI_iter_.fld';
% predfile = '\l2l2\TMVI_iter_.pre';

% mag_vecfile = 'Tile1_MVI.fld';
obsfile = 'MAG_Osborne_Sub.mag';

% mag_vecfile = '..\..\m_rem.dat';
% obsfile = '..\Obs_loc_TMI.obs';

% mag_vecfile = '..\magvec.fld';
% obsfile = '..\..\Obs_RAW_REM_GRID_TMI.obs';

mmax = 0.25;

cam_ang = [0 0];

vscale = 1;

% Define cutting planes
% nvec(1,1:3) = [0 -0.2 1]; xo = 1000; yo = 850  ; zo = 1420;
% nvec(2,1:3) = [0 1 0]; xo(2) = 1000; yo(2) = 730  ; zo(2) = 1440;
nvec(3,1:3) = [0 1 0]; x0 = 455970; y0 = 7556190  ; z0 = 300;

nvec = spdiags( 1./ sqrt(sum(nvec.^2,2)) , 0 ,3 ,3) * nvec;

cut_name = ['A','B','C'];

 
% Color scheme

%% Load in model and plot
set(figure, 'Position', [200 50 750 400]); 

ax1 = axes('Position',[0.55 .3 .425 .425]);

% Load geo unit
bif = load(geo_unit);
bif(:,3) = bif(:,3) - max(bif(:,3));

% Load data
[H, HI, HD, MI,MD, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir dsep obsfile]);
% H=1000;

[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

% Move to local coordinates
bif(:,1) = bif(:,1) - xn(1);
bif(:,2) = bif(:,2) - yn(1);
% bif(:,3) = bif(:,3) - zn(1);
% 
obsx = obsx - xn(1);
obsy = obsy - yn(1);
obsz = obsz - zn(1);

% x0 = x0 - xn(1);
% y0 = y0 - yn(1);
% z0 = z0 - zn(1);
% 
xn = xn - xn(1); 
yn = yn - yn(1); 
zn = zn - zn(1); 
% 
% 
% dx = xn(2:end) - xn(1:end-1); nx = length(dx);
% dy = yn(2:end) - yn(1:end-1); ny = length(dy);
% dz = zn(1:end-1) - zn(2:end); nz = length(dz);
% xx = (xn(2:end) + xn(1:end-1))/2;   xx = xx(padW+1:end-padE);
% yy = (yn(2:end) + yn(1:end-1))/2;   yy = yy(padS+1:end-padN);
% zz = (zn(2:end) + zn(1:end-1))/2;   zz = zz(padT+1:end-padB);
% 
% [XX,ZZ,YY] = meshgrid(xx,zz,yy); 

% Load mag true
% m_true = load([work_dir '\' model_true]);
% m_true = reshape(m_true,nz,nx,ny);


% Load magnetization vector

% ax1 = axes('Position',[0.15 -.1 0.7 0.7]);
%     h = imagesc(xx,yy,m2D); hold on

% scatter(min(XX2D(:)), max(YY2D(:)),'k.'); hold on
% text(min(XX2D(:)), max(YY2D(:)),['\textbf{' cut_name(ii) '}'],'interpreter', 'latex','FontSize',14,'VerticalAlignment','top');
% scatter(max(XX2D(:)), min(YY2D(:)),'k.')
% text(max(XX2D(:)), min(YY2D(:)),['\textbf{' cut_name(ii) '"}'],'interpreter', 'latex','FontSize',14,'HorizontalAlignment','right');
%     title('$Model$','interpreter', 'latex','FontSize',14);

%     qq = quiver(xx,yy,mx2D,my2D,'LineWidth',1,'Color','k','MaxHeadSize',1);


set(gca,'YDir','normal')

% for ii = 1 : 2
%     for jj = 1 : 2
%         j = (-1)^ii;
%         k = (-1)^jj;
%         quiver3(k*j*60,j*60,-140+j*60,0,-j*120,0,'LineWidth',2,'Color','k','ShowArrowHead','off','AutoScale','off');
%         quiver3(j*60,k*j*60,-140+j*60,-j*120,0,0,'LineWidth',2,'Color','k','ShowArrowHead','off','AutoScale','off');
%         quiver3(j*60,k*j*60,-140+j*60,0,0,-j*120,'LineWidth',2,'Color','k','ShowArrowHead','off','AutoScale','off');
%     end
% end


% quiver3(0,min(yy),-min(dz),0,-2*sum(yy),0,'LineWidth',2,'Color','k','ShowArrowHead','off','AutoScale','off');


% colormap('jet')
%plot_vec(ax1,XX2D,YY2D,ZZ2D,mx2D,my2D,mz2D,m2D,iso_cut_vec,0,0.5,1.25,3)

% axis([min(XX2D(:)) max(XX2D(:)) min(YY2D(:)) max(YY2D(:))])

% set(gca,'YTickLabel',[],'Box','on')
% set(gca,'XTickLabel',[],'Box','on')
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
% zlabh = get(gca,'ZLabel');
% set(zlabh,'Position',get(zlabh,'Position') + [50 0 0])

% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [1 0 0])

set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

% ax1.XTick = [x_ext(1)+25 x_ext(2)-25] ;
% ax1.XTickLabel = {num2str(x_ext(1)+25) ,num2str(x_ext(2)-25)};
% set(gca, 'YAxisLocation', 'right')
set(gca,'YDir','normal')
hold on
% title('$A-A"$','interpreter', 'latex','FontSize',14)
scatter3(bif(:,1),bif(:,2),bif(:,3),'k.','MarkerEdgeColor','k')

text(min(bif(:,1))+300,min(bif(:,2)),-320,'Ironstone','interpreter', 'latex','FontSize',12,'Color','k')
text(min(bif(:,1))+200,min(bif(:,2)),-30,'Mesozoic','interpreter', 'latex','FontSize',12,'Color','k')

quiver3(min(bif(:,1))+400,min(bif(:,2)),-300,100,0,100,'k','LineWidth',1,'MaxHeadSize',100)
quiver3(min(bif(:,1))+400,min(bif(:,2)),-300,-25,0,50,'k','LineWidth',1,'MaxHeadSize',100)

view(cam_ang)
cvec = mmax*[0 0.2 0.4 0.6 0.8 1.05];
bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax,'linear');
colormap(ax1,bb);
%colormap(ax2,bb);
% caxis(ax1,[0 mmax])
% ylim([y0-10, y0 + 10])
zlim([-350, 0])
xlim([min(bif(:,1)) max(bif(:,1))])
text(max(bif(:,1)),max(bif(:,2)),20, '$A"$','interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','Color','k')
text(min(bif(:,1)),min(bif(:,2)),20, '$A$','interpreter', 'latex','FontSize',14,'HorizontalAlignment','right','Color','k')

% caxis(ax2,[0 mmax])


%% Plot Data
% Load data
ax2 = axes('Position',[0.05 .25 .5 .5]);
nx = 100;
ny = 100;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / nx);
dy = (( ymax - ymin) / ny);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[X,Y] = meshgrid(x,y);

data_interp     = griddata(obsx, obsy, d,X,Y,'linear'); 
% data_interp(isnan(data_interp)) = min(data_interp(:))*2;


h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(ax2,jet);
contour(x,y,data_interp,'k'); hold on
plot([min(bif(:,1)) max(bif(:,1))],[min(bif(:,2)) max(bif(:,2))],'r--','LineWidth',2)
% scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
set(gca, 'XAxisLocation', 'top')
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% set(gca, 'YAxisLocation', 'right')
ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [100 0 0])
grid on
axis equal
axis tight

% quiver(1350,1500,150*(0.6428),150*(0.7660),'k','LineWidth',2,'MaxHeadSize',1);
% text(1350,1500,'$N$','interpreter', 'latex','FontSize',16,'Color','k','HorizontalAlignment','right','VerticalAlignment','top');
text(max(bif(:,1)),max(bif(:,2)), '$A"$','interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','Color','r')
text(min(bif(:,1)),min(bif(:,2)), '$A$','interpreter', 'latex','FontSize',14,'HorizontalAlignment','right','Color','r')
% title('$d^{Fwr}$', 'interpreter', 'latex','FontSize',14)
% text(min(xx)-dx*20, mean(yy),'$(a)$', 'interpreter', 'latex','FontSize',14)



%% Add colorbars
ax = axes('Position',[0.2 0.05 .20 .30]);
cbar = colorbar('SouthOutside');
colormap(ax,jet);
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([min(data_interp(:)) max(data_interp(:))]))
set(gca,'Visible','off');
text(0.5,-0.75,'$nT$', 'interpreter', 'latex','FontSize',12,'HorizontalAlignment','center')
text(0.4,-1.75,'$(a)$', 'interpreter', 'latex','FontSize',14)
text(2.75,-1.75,'$(b)$', 'interpreter', 'latex','FontSize',14)