% Create figure through 3D model
clear all
close all

%% INPUT PARAM
addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB;

ROT = [30 0 0];

nx = 3; ny = 3; nz = 3;

cam_ang =[30 45];
%%
xx = -1:1;
zz = -1:1;

dx = ones(nx,1);
dy = ones(ny,1);
dz = ones(nz,1);

nullcell = ones(nz*nx*ny,1);
mactv = sum(nullcell);

[A, GRAD, V] = get_GRAD_op3D_TENSIL_Kron(dx,dy,dz,nullcell,'FWR');

Rz = @(x)   [cosd(x) -sind(x) 0;
            sind(x) cosd(x) 0;
            0 0 1];

Ry = @(x)   [cosd(x) 0 -sind(x);
            0 1 0;
            sind(x) 0 cosd(x)];

Rx = @(x)   [1 0 0;
            0 cosd(x) -sind(x);
            0 sind(x) cosd(x)];
        
rz = Rz(ROT(1));
%     rz = rz*spdiags(sum(abs(rz),2).^-1,0,3,3);

ry = Ry(ROT(2));
%     ry = ry*spdiags(sum(abs(ry),2).^-1,0,3,3);

rx = Rx(ROT(3));

rot = rx*ry*rz;
% Scale the rows


% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(1,:)'),'descend') ;

Gx = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while denom < (45 - 1e-4)
    
    Gx = Gx + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gx = Gx * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gx),2)) ~= 2;
% Gx(indx,:) = 0;

% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(2,:)'),'descend') ;

Gy = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while denom < (45 - 1e-4)
    
    Gy = Gy + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gy = Gy * spdiags(ones(mactv,1)/denom,0,mactv,mactv);
% indx = round(sum(abs(Gy),2)) ~= 2;
% Gy(indx,:) = 0;


% Get index and weight for gradients Gx
[val,ind] = sort(45 - acosd(A*rot(3,:)'),'descend') ;

Gz = GRAD{ind(1)} * val(1);
denom = val(1);
count = 2;
while denom < (45 - 1e-4)
    
    Gz = Gz + GRAD{ind(count)} * val(count);
    denom = denom + val(count);
    count = count + 1;
    
end

Gz = Gz * spdiags(ones(mactv,1)/denom,0,mactv,mactv);

%% Create model and plot gradients
m = zeros(nz,nx,ny);
m(2,2,2) = 1;
m = m(:);

gradx = Gx * m;
grady = Gy * m;
gradz = Gz * m;

grad = gradx + grady;
m2D = reshape(grad,nz,nx,ny);
m2D = squeeze( m2D(2 , : , :));
[XX2D,ZZ2D] = ndgrid(xx,zz);
YY2D = ones(size(XX2D))*1;
h = imagesc(xx,zz,m2D'); hold on

% m2D = reshape(gradx,nz,nx,ny);
% m2D = squeeze( m2D(: , : , 2));
% [XX2D,ZZ2D] = ndgrid(xx,zz);
% YY2D = ones(size(XX2D))*0;
% h = surf(YY2D,ZZ2D,XX2D,m2D); hold on
% 
% m2D = reshape(gradx,nz,nx,ny);
% m2D = squeeze( m2D(: , : , 1));
% [XX2D,ZZ2D] = ndgrid(xx,zz);
% YY2D = ones(size(XX2D))*-1;
% h = surf(YY2D,ZZ2D,XX2D,m2D); hold on

alpha(1)
view(cam_ang)
set(gca,'YDir','normal')

% axis([-1 1 -1 1 -1 1])
set(gca,'ZTickLabel',[],'Box','on')
set(gca,'YTickLabel',[],'Box','on')
set(gca,'XTickLabel',[],'Box','on')
zlabel('$z$', 'interpreter', 'latex','FontSize',14)
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
zlabh = get(gca,'ZLabel');

ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') - [0.5 0 0])

xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [0 0.5 0])

set(get(gca,'YLabel'),'Rotation',360);
set(get(gca,'ZLabel'),'Rotation',360);
axis equal
grid on

% set(gca, 'YAxisLocation', 'right')
set(gca,'YDir','normal')
hold on


cvec = [-1.1 0.0 1.1];
bb = interp1([cvec'],[0 100 255;255 255 255;255 0 0]/255,-1.1:1e-4:1.1,'linear');
colormap(bb);
% colorbar;
%colormap(ax2,bb);
caxis([-1 1])
% caxis(ax2,[0 mmax])

%% ADD VECTORS
qq = quiver(0,0,rot(1,1),rot(1,2),'Marker','o','MaxHeadSize',0.5,'LineWidth',2,'Color','k');
qq = quiver(0,0,rot(2,1),rot(2,2),'Marker','o','MaxHeadSize',0.5,'LineWidth',2,'Color','k');
qq = quiver3(0,0,0,0,0,1,'Marker','o','MaxHeadSize',0.5,'LineWidth',2,'Color','k');
%% Add color bar

plot([0 1],[0 1],':k','LineWidth',1)
plot([0 1],[0 -1],':k','LineWidth',1)
plot([0 0],[0 1],':k','LineWidth',1)
plot([0 1],[0 0],':k','LineWidth',1)
plot([-1 0],[0 0],':k','LineWidth',1)
plot([0 -1],[0 1],':k','LineWidth',1)
plot([0 -1],[0 -1],':k','LineWidth',1)
plot([0 -1],[0 -1],':k','LineWidth',1)
plot([0 0],[0 -1],':k','LineWidth',1)
plot3([0 0],[0 1],[0 1],'k:')
plot3([0 0],[0 0],[0 -1],'k:')
plot3([0 -1],[0 -1],[0 -1],'k:')
plot3([0 1],[0 0],[0 -1],'k:')
plot3([0 -1],[0 0],[0 1],'k:')
plot3([0 -1],[0 -1],[0 1],'k:')
plot3([0 1],[0 -1],[0 -1],'k:')
plot3([0 0],[0 -1],[0 -1],'k:')
plot3([0 1],[0 1],[0 -1],'k:')

plot([-.5 -.5],[-.5 1.5],'k','LineWidth',2)
plot([-.5 1.5],[1.5 1.5],'k','LineWidth',2)
plot([1.5 1.5],[1.5 0.5],'k','LineWidth',2)
plot([1.5 0.5],[0.5 0.5],'k','LineWidth',2)
plot([0.5 0.5],[0.5 -0.5],'k','LineWidth',2)
plot([-0.5 0.5],[-0.5 -0.5],'k','LineWidth',2)

plot([-.5 -.5],[-.5 0.5],'k--','LineWidth',2)
plot([-.5 1.5],[0.5 0.5],'k--','LineWidth',2)
plot([1.5 1.5],[0.5 -1.5],'k--','LineWidth',2)
plot([1.5 0.5],[-1.5 -1.5],'k--','LineWidth',2)
plot([0.5 0.5],[-1.5 -0.5],'k--','LineWidth',2)
plot([-0.5 0.5],[-0.5 -0.5],'k--','LineWidth',2)

text(1,0,0.1,'$g_x$', 'interpreter', 'latex','FontSize',12)
text(0,1,0.1,'$g_y$', 'interpreter', 'latex','FontSize',12)
text(1,1,0.1,'$g_{yx}$', 'interpreter', 'latex','FontSize',12)
text(1,-1,0.1,'$g_{xy}$', 'interpreter', 'latex','FontSize',12)
text(-.1,0,1.05,'$\hat G_{z}$', 'interpreter', 'latex','FontSize',14)
text(rot(1,1),rot(1,2),0.1,'$\hat G_{y}$', 'interpreter', 'latex','FontSize',14)
text(rot(2,1),rot(2,2),0.1,'$\hat G_{x}$', 'interpreter', 'latex','FontSize',14)

axis([-1.5 1.5 -1.5 1.5 -1 1]);
axis square