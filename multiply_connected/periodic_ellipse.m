function obstacle = periodic_ellipse(centre,semiaxes,angle)
% Elliptical solid; semiaxes=[a,b], rotation in radians, physical coordinates.
if nargin<3, angle=0; end
validateattributes(centre,{'numeric'},{'scalar','finite'});
validateattributes(semiaxes,{'numeric'},{'real','positive','finite','numel',2});
validateattributes(angle,{'numeric'},{'real','scalar','finite'});
a=semiaxes(1); b=semiaxes(2); rotation=exp(1i*angle);
obstacle=struct('type','ellipse','centre',centre,'semiaxes',[a,b],'angle',angle, ...
    'marker',centre,'boundary',@(t) centre+rotation*(a*cos(t)+1i*b*sin(t)), ...
    'derivative',@(t) rotation*(-a*sin(t)+1i*b*cos(t)));
end
