function cuts = periodic_obstacle_sections(o,x)
% Solid vertical intervals; ellipses use an exact quadratic intersection.
cuts=zeros(0,2); dx=mod(x-real(o.marker)+pi,2*pi)-pi;
if strcmp(o.type,'ellipse')
    dx=real(o.marker)+dx-real(o.centre);
    a=o.semiaxes(1); b=o.semiaxes(2); c=cos(o.angle); s=sin(o.angle);
    A=s^2/a^2+c^2/b^2; B=2*dx*c*s*(1/a^2-1/b^2);
    C=dx^2*(c^2/a^2+s^2/b^2)-1; discriminant=B^2-4*A*C;
    if discriminant>0
        cuts=imag(o.centre)+(-B+[-1,1]*sqrt(discriminant))/(2*A);
    end
    return
end
% Find all resolved transverse crossings, then refine on the original curve.
% Vertical tangencies contribute no solid interval; very small unresolved
% lobes require denser geometry sampling, as with polygon membership.
xx=real(o.marker)+dx; m=numel(o.polygon); t=2*pi*(0:m)'/m;
values=real(o.boundary(t))-xx; roots=[];
for k=1:m
    if values(k)==0, roots(end+1)=t(k); end %#ok<AGROW>
    if values(k)*values(k+1)<0
        roots(end+1)=fzero(@(s) real(o.boundary(s))-xx,[t(k),t(k+1)]); %#ok<AGROW>
    end
end
if isempty(roots), return; end
ys=sort(imag(o.boundary(roots))); ys=ys([true,diff(ys)>1e-11]);
for k=1:numel(ys)-1
    if periodic_obstacle_contains(xx+1i*mean(ys(k:k+1)),o)
        cuts(end+1,:)=ys(k:k+1); %#ok<AGROW>
    end
end
end
