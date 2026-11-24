function [ax,data] = plot_periodic_stokes_prsa(r,ax,data)
% PRSA-style combined speed, streamline and pole plot over two periods.
% All AAA poles are red; Laurent-series centres are blue dots (RSOS Fig. 8).
% The display window matches poiseuille_smooth.m for its channel geometries.
% All retained poles remain in the solver, including poles outside the view.
% Optional data contains one-period x,y,speed,psi grids (psi=0 on bottom wall).
if nargin<2 || isempty(ax),figure('Color','w');ax=axes;end
g=r.geometry;
if nargin<3 || isempty(data)
 x=linspace(0,2*pi,401);[top,bot]=periodic_channel_walls(g,x);
 y=linspace(min(imag(bot)),max(imag(top)),161);
 [xx,yy]=meshgrid(x,y);z=xx+1i*yy;inside=r.inFluid(z);
 speed=nan(size(z));ps=speed;ix=find(inside);
 [~,bottom]=periodic_channel_walls(g,0);p0=r.psi(bottom);
 for first=1:2000:numel(ix)
  ii=ix(first:min(first+1999,numel(ix)));
  speed(ii)=abs(r.uv(z(ii)));ps(ii)=r.psi(z(ii))-p0;
 end
 data=struct('x',x,'y',y,'speed',speed,'psi',ps);
end
% Repeat periodic fields exactly, avoiding duplicate samples at the seams.
x=[data.x-2*pi,data.x(2:end),data.x(2:end)+2*pi];
speed=[data.speed,data.speed(:,2:end),data.speed(:,2:end)];
ps=[data.psi,data.psi(:,2:end),data.psi(:,2:end)];
keep=x>=-pi & x<=3*pi;x=x(keep);speed=speed(:,keep);ps=ps(:,keep);
y=data.y;[top,bot]=periodic_channel_walls(g,x);
hold(ax,'on');pcolor(ax,x,y,speed);shading(ax,'interp');colormap(ax,parula);
caxis(ax,[0 max(speed(:),[],'omitnan')]);
cb=colorbar(ax);cb.Label.String='Velocity magnitude';cb.FontSize=26;cb.Label.FontSize=26;
Q=periodic_holes_flux(r,0);contour(ax,x,y,ps,sort(Q*(.1:.1:.9)),'k','LineWidth',.6);
pmin=min(ps(:),[],'omitnan');
if pmin < -1e-9
 contour(ax,x,y,ps,sort(pmin*[.02 .06 .15 .3 .5 .7 .9]),'y','LineWidth',.7);
end
% Upper-wall eddies have streamfunction values above the throughflow range.
pmax=max(ps(:),[],'omitnan');
if Q>0 && pmax>Q+1e-9
 contour(ax,x,y,ps,Q+(pmax-Q)*[.02 .06 .15 .3 .5 .7 .9],'y','LineWidth',.7);
end
% Optional diagnostics resolve tiny counter-rotating bottom eddies in Couette flow.
if isfield(data,'criticalPoints')
 critical=data.criticalPoints;
 for j=1:size(critical.xy,1)
  pc=critical.psi(j);yc=critical.xy(j,2);
  if pc>1e-9 && pc<abs(Q)*1e-3 && critical.determinant(j)>0
   local=ps;local(y>yc+.4,:)=NaN;
   contour(ax,x,y,local,pc*[.15 .35 .55 .75 .9],'w','LineWidth',.8);
  end
 end
end
plot(ax,real(top),imag(top),'k','LineWidth',.8);
plot(ax,real(bot),imag(bot),'k','LineWidth',.8);
t=linspace(0,2*pi,501);
for copy=-1:1
 for j=1:numel(g.centres)
  o=periodic_obstacle_at(g,j);curve=o.boundary(t)+copy*2*pi;
  fill(ax,real(curve),imag(curve),[.88 .88 .88],'EdgeColor','k','LineWidth',.8);
  if ~isempty(r.basis.holeAAA{j})
   p=r.basis.holeAAA{j}.kept+copy*2*pi;
   plot(ax,real(p),imag(p),'.r','MarkerSize',8);
  end
  plot(ax,real(o.marker)+copy*2*pi,imag(o.marker),'.b','MarkerSize',12);
 end
 for j=1:numel(r.basis.wallAAA)
  p=r.basis.wallAAA(j).kept+copy*2*pi;
  plot(ax,real(p),imag(p),'.r','MarkerSize',8);
 end
end
% Match the original PRSA Fig. 3 vertical extent; enlarge for taller channels.
yl=[min(-1.6,min(imag(bot))-.3),max(1.6,max(imag(top))+.3)];
plot(ax,[0 0],yl,'k--','LineWidth',1.2);
plot(ax,[2*pi 2*pi],yl,'k--','LineWidth',1.2);
axis(ax,'equal');axis(ax,'off');xlim(ax,[-pi 3*pi]);ylim(ax,yl);
end
