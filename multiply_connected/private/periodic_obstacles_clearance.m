function [gap,clearance] = periodic_obstacles_clearance(g,width)
% Smooth-curve numerical geometry checks, independent of the collocation grid.
% Wall clearances are vertical; inter-obstacle gaps are sampled/refined distances.
J=numel(g.obstacles); gap=width; clearance=zeros(J,2); m=4096;
t=2*pi*(0:m-1)'/m; h=2*pi/m;
for j=1:J
    o=g.obstacles{j}; z=o.boundary(t);
    for side=1:2
        f=@(theta) wall_gap(g,o,theta,side); v=f(t); minimum=min(v);
        ix=find(v<=circshift(v,1) & v<=circshift(v,-1) ...
            & (v<circshift(v,1) | v<circshift(v,-1)));
        for k=ix(:).'
            [~,d]=fminbnd(f,t(k)-h,t(k)+h,optimset('TolX',1e-12,'Display','off'));
            minimum=min(minimum,d);
        end
        clearance(j,side)=minimum;
    end
    gap=min([gap,clearance(j,:),2*pi-(max(real(z))-min(real(z)))]);
    assert(gap>0,'Obstacle touches/intersects a wall or its periodic copy.');
    for k=1:j-1
        other=g.obstacles{k}; w=other.boundary(t);
        assert(~any(periodic_obstacle_contains(z,other)) && ~any(periodic_obstacle_contains(w,o)), ...
            'Obstacles overlap or contain one another.');
        % Coarse nearest pair, then continuous minimization of both parameters.
        ts=t(1:4:end); zs=o.boundary(ts); ws=other.boundary(ts);
        dx=mod(real(zs)-real(ws.')+pi,2*pi)-pi;
        distance=dx.^2+(imag(zs)-imag(ws.')).^2;
        [~,idx]=min(distance(:)); [a,b]=ind2sub(size(distance),idx);
        [~,d]=fminsearch(@(ab) pair_distance(o,other,ab),[ts(a),ts(b)], ...
            optimset('TolX',1e-12,'TolFun',1e-24,'Display','off'));
        assert(d>1e-20,'Obstacle boundaries touch or intersect.'); gap=min(gap,sqrt(d));
    end
end
end
function d = wall_gap(g,o,t,side)
z=o.boundary(t); [top,bot]=periodic_channel_walls(g,real(z));
if side==1, d=imag(top-z); else, d=imag(z-bot); end
end
function d = pair_distance(a,b,t)
z=a.boundary(t(1))-b.boundary(t(2)); dx=mod(real(z)+pi,2*pi)-pi;
d=dx^2+imag(z)^2;
end
