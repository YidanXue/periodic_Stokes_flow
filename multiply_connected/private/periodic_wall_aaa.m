function [Pol,metadata] = periodic_wall_aaa(g,top,bot,tol,mmax)
% The original Schwarz-function AAA construction, separately on each wall.
% Physical poles are -i*log(zeta_pole), wrapped to [0,2*pi).
% Replace the artificial-cell polygon filter by exact graph/periodic disk
% membership. Only poles strictly in excluded material are admitted.
assert(exist('aaa','file')~=0, ...
    'Curved-wall AAA requires Chebfun aaa.m on the MATLAB path.');
Pol = {}; metadata = struct('source',{},'rawMapped',{},'kept',{}, ...
    'rejectedFluidOrBoundary',{},'discardedNonfinite',{},'errorHistory',{});
for side = 1:2
    if side==1, wall=top; name='top'; else, wall=bot; name='bottom'; end
    zeta = exp(1i*wall);
    [~,p,~,~,~,~,~,err] = aaa(conj(zeta),zeta,'tol',tol,'mmax',mmax);
    valid = isfinite(p) & p~=0; removed = sum(~valid); p = -1i*log(p(valid));
    % Preserve original ordering: nonnegative real parts, then negative ones.
    p = [p(real(p)>=0);p(real(p)<0)+2*pi];
    p = mod(real(p),2*pi)+1i*imag(p);
    [t,b] = periodic_channel_walls(g,real(p));
    guard = 100*eps*(1+abs(imag(p)));
    outside = imag(p)>imag(t)+guard | imag(p)<imag(b)-guard;
    for j = 1:numel(g.centres)
        outside = outside | periodic_obstacle_contains(p,periodic_obstacle_at(g,j),true);
    end
    kept = p(outside);
    if ~isempty(kept), Pol{end+1} = kept; end %#ok<AGROW>
    metadata(side) = struct('source',name,'rawMapped',p,'kept',kept, ...
        'rejectedFluidOrBoundary',p(~outside),'discardedNonfinite',removed,'errorHistory',err);
end
end
