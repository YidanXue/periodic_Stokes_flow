function [Z,component,weights] = periodic_holes_boundary(geom,mWall,mHole,shift)
% Half-open, optionally shifted samples; obstacle coordinates may cross the seam.
if nargin < 4, shift = 0; end
J = numel(geom.centres);
if isscalar(mHole), mHole = repmat(mHole,1,J); end
x = 2*pi*((0:mWall-1)'+shift)/mWall;
[top,bot] = periodic_channel_walls(geom,x);
Z = [top; flipud(bot)];
component = [ones(mWall,1); 2*ones(mWall,1)];
weights = repmat(2*pi/mWall,2*mWall,1);
if isfield(geom,'h1')
    % Arclength weights are metadata; the fit remains unweighted as originally.
    h = 1e-5;
    [tp,bp] = periodic_channel_walls(geom,x+h);
    [tm,bm] = periodic_channel_walls(geom,x-h);
    weights = (2*pi/mWall)*[abs((tp-tm)/(2*h));flipud(abs((bp-bm)/(2*h)))];
end
for j = 1:J
    theta = 2*pi*((0:mHole(j)-1)'+shift)/mHole(j);
    o=periodic_obstacle_at(geom,j);
    Z = [Z; o.boundary(theta)]; %#ok<AGROW>
    component = [component; (j+2)*ones(mHole(j),1)]; %#ok<AGROW>
    weights = [weights; (2*pi/mHole(j))*abs(o.derivative(theta))]; %#ok<AGROW>
end
end
