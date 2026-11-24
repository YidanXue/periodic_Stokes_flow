function o = periodic_obstacle_at(g,j)
% Common access without changing the original centres/radii interface.
if isfield(g,'obstacles'), o=g.obstacles{j};
else, o=periodic_ellipse(g.centres(j),[g.radii(j),g.radii(j)]); end
end
