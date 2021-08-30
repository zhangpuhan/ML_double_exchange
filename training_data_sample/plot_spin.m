
clear all;
close all;
format long;
format compact;

color=dlmread('rainbow.dat');

% Read data
m=dlmread('c15800.dat');


x=m(:,1);
y=m(:,2);
Sx=m(:,4);
Sy=m(:,5);
Sz=m(:,6);
hx=m(:,7);
hy=m(:,8);
hz=m(:,9);
nd=m(:,3);
bx=m(:,3);

xr=x;
yr=y;
for i=1:length(x)
   xr(i) = x(i) - 0.5 * Sx(i);
   yr(i) = y(i) - 0.5 * Sz(i);
end

ix=find(diff(x)<0,1,'first');
iy=length(x)/ix;
% FLIP=repmat((1:ix)',1,iy)+repmat((1:iy),ix,1);
% FLIP=1-2*mod(FLIP,2);
% FLIP=reshape(FLIP,[],1);
% dx=dx.*FLIP;
% dy=dy.*FLIP;

% Vector -> matrix
xm=reshape(x,ix,iy);
ym=reshape(y,ix,iy);
zm=reshape(bx,ix,iy);

F = griddedInterpolant(xm,ym,zm,'linear','cubic');
[xq,yq] = ndgrid(0:0.05:29,0:0.05:29);
vq = F(xq, yq);


% % % Plot figure
% fi=figure('Position',[200,150,734,550]);
% ax=axes('Position',[0.12,0.12,0.80,0.80]);
% hold(ax,'on');
% % % %surf(ax,xm,ym,zm);
% % % %pcolor(ax,xm,ym,zm);
% % % 
% surf(ax,xq, yq, vq);
% pcolor(ax,xq,yq,vq);
% shading(ax,'interp')
% caxis([0.2 1])
% % % 
% quiver(ax,xr,yr,Sx,Sz,'Color','white','AutoScaleFactor',0.45,'LineWidth',0.9,'MaxHeadSize',0.85);
% % % 
% hold(ax,'off');
% shading(ax,'interp');
% % % 
% % % % colormap(ax,'jet');
% colormap(gca,'jet');
% colorbar;
% set(gca,'Box','on','FontSize',20,'TickLabelInterpreter','latex')
% % % %set(gca,'XLim',[min(x),max(x)],'YLim',[min(y),max(y)]);
% set(gca,'XTick',0:5:59,'YTick',0:5:59);
% axis equal;
% % % 
% % % %xlabel('$x$','FontSize',20,'Interpreter','latex');
% % % %ylabel('$y$','FontSize',20,'Interpreter','latex');
% % % 
% saveas(fi,'spin.png');


z2m=reshape(nd,ix,iy);

G = griddedInterpolant(xm,ym,z2m,'cubic','nearest');
[xq,yq] = ndgrid(0:0.05:29,0:0.05:29);
vq = G(xq, yq);

% Plot figure
fi2=figure('Position',[150,150,700,640]);
ax2=axes('Position',[0.12,0.12,0.80,0.80]);
colorbar;

hold(ax2,'all');
surf(ax2,xq,yq,vq);
%surf(ax2,xm,ym,z2m);
%pcolor(ax2,xm,ym,z2m);
hold(ax2,'off');
shading(ax2,'interp');
caxis([0.35 0.5])

colormap(ax2,'jet');
set(gca,'Box','on','FontSize',20,'TickLabelInterpreter','latex')
%set(gca,'XLim',[min(x),max(x)],'YLim',[min(y),max(y)]);
set(gca,'XTick',0:5:29,'YTick',0:5:29);
axis equal;

%xlabel('$x$','FontSize',20,'Interpreter','latex');
%ylabel('$y$','FontSize',20,'Interpreter','latex');

saveas(fi2,'nd_15800.png');


