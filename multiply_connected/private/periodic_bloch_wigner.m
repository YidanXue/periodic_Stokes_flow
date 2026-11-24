function D = periodic_bloch_wigner(w)
% Real, single-valued Bloch-Wigner dilogarithm; no symbolic toolbox.
% D(w)=-D(1/w)=-D(1-w). Reduce to |w|<=1, Re(w)<=1/2,
% then use Li2(1-exp(-t)) = t-t^2/4+sum B_2k*t^(2k+1)/(2k+1)!.
% In the reduced region |t| is far below the 2*pi convergence radius.
sz = size(w); w = w(:); signs = ones(size(w));
ix = abs(w)>1; w(ix) = 1./w(ix); signs(ix) = -signs(ix);
ix = real(w)>0.5; w(ix) = 1-w(ix); signs(ix) = -signs(ix);
D = zeros(size(w)); ix = w~=0;
z = w(ix); t = -log(1-z); L = t-t.^2/4;
B = [1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510,43867/798,-174611/330];
for k = 1:numel(B), L = L+B(k)*t.^(2*k+1)/factorial(2*k+1); end
D(ix) = signs(ix).*(imag(L)+log(abs(z)).*angle(1-z));
D = reshape(D,sz);
end
