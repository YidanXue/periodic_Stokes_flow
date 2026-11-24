function [Q,intervals,pieces] = periodic_holes_flux(result,x,nQuad)
% Integrate u over fluid intervals only, using Gauss-Legendre quadrature.
if nargin < 3, nQuad = 160; end
validateattributes(x,{'numeric'},{'real','scalar','finite'});
validateattributes(nQuad,{'numeric'},{'real','scalar','integer','>=',2});
g = result.geometry; cuts = zeros(0,2);
[top,bot] = periodic_channel_walls(g,x);
for j = 1:numel(g.centres)
    cuts=[cuts;periodic_obstacle_sections(periodic_obstacle_at(g,j),x)]; %#ok<AGROW>
end
cuts = sortrows(cuts,1); intervals = zeros(0,2); lo = imag(bot);
for j = 1:size(cuts,1)
    if cuts(j,1)>lo, intervals(end+1,:) = [lo,cuts(j,1)]; end %#ok<AGROW>
    lo = max(lo,cuts(j,2));
end
if lo<imag(top), intervals(end+1,:) = [lo,imag(top)]; end
k = (1:nQuad-1)'; beta = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(beta,1)+diag(beta,-1));
[nodes,ix] = sort(diag(D)); weights = 2*V(1,ix)'.^2;
pieces = zeros(size(intervals,1),1);
for j = 1:numel(pieces)
    mid = mean(intervals(j,:)); half = diff(intervals(j,:))/2;
    pieces(j) = half*sum(weights.*real(result.uv(x+1i*(mid+half*nodes))));
end
Q = sum(pieces);
end
