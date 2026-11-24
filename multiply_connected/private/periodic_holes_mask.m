function inside = periodic_holes_mask(z,geom)
% Fluid membership, including all periodic copies of every solid.
[top,bot] = periodic_channel_walls(geom,real(z));
inside = imag(z)>imag(bot) & imag(z)<imag(top);
for j = 1:numel(geom.centres)
    inside = inside & ~periodic_obstacle_contains(z,periodic_obstacle_at(geom,j));
end
end
