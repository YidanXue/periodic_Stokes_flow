function g = periodic_obstacles_prepare(g)
% Optional geom.obstacles is a cell array of smooth simple closed boundaries.
% Each custom obstacle needs boundary(t), marker, and optionally derivative(t).
if ~isfield(g,'obstacles'), return; end
if isstruct(g.obstacles), g.obstacles=num2cell(g.obstacles); end
assert(iscell(g.obstacles),'obstacles must be a cell or struct array.');
assert(isempty(g.radii),'Use either obstacles or centres/radii, not both.');
g.centres=zeros(1,numel(g.obstacles));
t=2*pi*(0:4095)'/4096;
for j=1:numel(g.obstacles)
    o=g.obstacles{j};
    assert(isfield(o,'boundary') && isa(o.boundary,'function_handle') && isfield(o,'marker'), ...
        'Each obstacle needs a boundary(t) and an interior marker.');
    validateattributes(o.marker,{'numeric'},{'scalar','finite'});
    if ~isfield(o,'type'), o.type='parametric'; end
    z=o.boundary(t);
    assert(isequal(size(z),size(t)) && all(isfinite(z)), 'Boundary must be finite and vectorized.');
    assert(max(abs(o.boundary(t+2*pi)-z))<1e-10,'Obstacle boundary must be 2*pi-periodic in its parameter.');
    assert(max(abs(real(z-o.marker)))<pi,'Obstacle must fit inside the marker-centred horizontal period.');
    if ~isfield(o,'derivative')
        boundary=o.boundary;
        o.derivative=@(t) (boundary(t+1e-5)-boundary(t-1e-5))/(2e-5);
    end
    dz=o.derivative(t);
    assert(all(isfinite(dz)) && min(abs(dz))>1e-10,'Boundary parameterization must be regular.');
    o.polygon=z; o.bounds=[min(real(z)),max(real(z)),min(imag(z)),max(imag(z))];
    o.maskTolerance=2*max(abs(o.boundary(t+pi/4096)-(z+circshift(z,-1))/2))+1e-12;
    [in,on]=inpolygon(real(o.marker),imag(o.marker),real(z),imag(z));
    assert(in && ~on,'Completion marker must be strictly inside its obstacle.');
    if ~strcmp(o.type,'ellipse')
        assert(simple_polygon(z),'Obstacle boundary must be simple, without self-intersections.');
    end
    g.obstacles{j}=o; g.centres(j)=o.marker;
end
end

function simple = simple_polygon(z)
% Detect crossings of nonadjacent polygon edges, including collinear contacts.
n=numel(z); a=z; b=circshift(z,-1); simple=true;
for i=1:n-2
    js=(i+2:n)'; if i==1, js(js==n)=[]; end
    aa=a(js); bb=b(js);
    cross=@(u,v) imag(conj(u).*v);
    s1=cross(b(i)-a(i),aa-a(i)); s2=cross(b(i)-a(i),bb-a(i));
    t1=cross(bb-aa,a(i)-aa); t2=cross(bb-aa,b(i)-aa);
    boxes=max(min(real(a(i)),real(b(i))),min(real(aa),real(bb))) ...
        <=min(max(real(a(i)),real(b(i))),max(real(aa),real(bb))) ...
        & max(min(imag(a(i)),imag(b(i))),min(imag(aa),imag(bb))) ...
        <=min(max(imag(a(i)),imag(b(i))),max(imag(aa),imag(bb)));
    if any(boxes & s1.*s2<=0 & t1.*t2<=0), simple=false; return; end
end
end
