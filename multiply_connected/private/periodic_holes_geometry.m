function info = periodic_holes_geometry(g)
% Check smooth periodic graph walls and disjoint circular or parametric solids.
% Clearances are numerical estimates: sample 4096 points, then refine every
% sampled local minimum. Arbitrarily unresolved wall oscillations are excluded.
assert(isscalar(g.P) && isfinite(g.P) && abs(g.P-2*pi)<1e-14,'Only P=2*pi is implemented.');
assert(isreal(g.walls) && numel(g.walls)==2 && all(isfinite(g.walls)) && diff(g.walls)>0,'Invalid walls.');
if ~isfield(g,'obstacles')
assert(numel(g.centres)==numel(g.radii) && all(isfinite(g.centres)) ...
    && isreal(g.radii) && all(isfinite(g.radii) & g.radii>0),'Invalid circles.');
end
curved = isfield(g,'h1') || isfield(g,'h2');
assert(~curved || (isfield(g,'h1') && isfield(g,'h2')), 'Supply both h1 and h2.');
m = 4096; x = 2*pi*(0:m-1)'/m;
if curved
    for fun = {g.h1,g.h2}
        f = fun{1}; assert(isa(f,'function_handle'),'h1 and h2 must be function handles.');
        v = f(x); assert(isequal(size(v),size(x)) && all(isfinite(v)), 'Walls must be finite, vectorized graph functions.');
        assert(max(abs(real(v)-x))<1e-11,'Wall parameter must equal its horizontal coordinate.');
        for shift = [-2*pi,2*pi]
            assert(max(abs(f(x+shift)-v-shift))<1e-10, 'Wall heights must be 2*pi periodic.');
        end
    end
end
[top,bot] = periodic_channel_walls(g,x);
if curved
    width = periodic_min(@(s) height_gap(g,s),x);
else
    width = diff(g.walls);
end
assert(width>0,'Top and bottom walls must not touch or intersect.');
gap = width; wallClearance = zeros(numel(g.centres),2);
if isfield(g,'obstacles')
    [gap,wallClearance]=periodic_obstacles_clearance(g,width);
else
for j = 1:numel(g.centres)
    c = g.centres(j); r = g.radii(j); [t,b] = periodic_channel_walls(g,real(c));
    assert(imag(c)>imag(b) && imag(c)<imag(t),'A circle centre lies outside the channel.');
    if curved
        for side = 1:2
            wallClearance(j,side) = sqrt(periodic_min(@(s) distance_squared(g,c,s,side),x))-r;
        end
    else
        wallClearance(j,:) = [g.walls(2)-imag(c)-r,imag(c)-g.walls(1)-r];
    end
    gap = min([gap,wallClearance(j,:),2*pi-2*r]);
    for k = 1:j-1
        dx = mod(real(c-g.centres(k))+pi,2*pi)-pi;
        gap = min(gap,hypot(dx,imag(c-g.centres(k)))-r-g.radii(k));
    end
end
end
assert(gap>0,'Circles must be disjoint from curved walls and all periodic copies.');
info = struct('minimumGap',gap,'wallClearance',wallClearance, ...
    'minimumVerticalWidth',width,'sampleCount',m,'curved',curved, ...
    'heightRange',[min(imag(bot)),max(imag(top))]);
end

function value = height_gap(g,x)
[t,b] = periodic_channel_walls(g,x); value = imag(t-b);
end
function value = distance_squared(g,c,x,side)
[t,b] = periodic_channel_walls(g,x); if side==1, wall=t; else, wall=b; end
dx = mod(x-real(c)+pi,2*pi)-pi;
value = dx.^2+(imag(wall)-imag(c)).^2;
end
function value = periodic_min(f,x)
v = f(x); value = min(v); h = 2*pi/numel(x);
% Flat plateaus need no refinement; otherwise refine all resolved local minima.
ix = find(v<=circshift(v,1) & v<=circshift(v,-1) ...
    & (v<circshift(v,1) | v<circshift(v,-1)));
for k = ix(:).'
    [~,local] = fminbnd(f,x(k)-h,x(k)+h,optimset('TolX',1e-12,'Display','off'));
    value = min(value,local);
end
end
