function inside = periodic_obstacle_contains(z,o,strict)
% Periodic membership. Ellipses are analytic; custom curves use a cached polygon.
if nargin<3, strict=false; end
local=mod(real(z-o.marker)+pi,2*pi)-pi+1i*imag(z-o.marker);
if strcmp(o.type,'ellipse')
    w=(o.marker+local-o.centre)*exp(-1i*o.angle);
    value=(real(w)/o.semiaxes(1)).^2+(imag(w)/o.semiaxes(2)).^2;
    if strict, inside=value<1-1e-10; else, inside=value<=1; end
else
    zz=o.marker+local; inside=false(size(z)); box=o.bounds;
    ix=real(zz)>=box(1) & real(zz)<=box(2) & imag(zz)>=box(3) & imag(zz)<=box(4);
    [in,on]=inpolygon(real(zz(ix)),imag(zz(ix)),real(o.polygon),imag(o.polygon));
    if strict
        in=in & ~on;
        % Keep poles farther inside than the estimated polygon chord error.
        queries=zz(ix); a=o.polygon; v=circshift(a,-1)-a;
        for k=find(in(:)).'
            s=max(0,min(1,real((queries(k)-a).*conj(v))./abs(v).^2));
            if min(abs(queries(k)-a-s.*v))<=o.maskTolerance, in(k)=false; end
        end
    end
    inside(ix)=in;
end
end
