function [p,metadata] = periodic_obstacle_aaa(o,Z,tol,mmax)
% Original Schwarz AAA in zeta=exp(i*z), retaining poles in THIS solid.
assert(exist('aaa','file')~=0,'Obstacle AAA requires Chebfun aaa.m on the path.');
zeta=exp(1i*Z);
[~,raw,~,~,~,~,~,err]=aaa(conj(zeta),zeta,'tol',tol,'mmax',mmax);
valid=isfinite(raw) & raw~=0; discarded=sum(~valid); raw=-1i*log(raw(valid));
raw=real(o.marker)+mod(real(raw)-real(o.marker)+pi,2*pi)-pi+1i*imag(raw);
keep=periodic_obstacle_contains(raw,o,true); p=raw(keep);
metadata=struct('rawMapped',raw,'kept',p,'rejected',raw(~keep), ...
    'discardedNonfinite',discarded,'errorHistory',err);
end
