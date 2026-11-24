function [W,P,O] = periodic_hole_columns(z,centres)
% Branch-free completion; W=u-i*v, p-i*omega=4*f'. Period 2*pi, mu=1.
% Columns for each solid are REAL coefficients [alpha,beta,C].
z = z(:); W = zeros(numel(z),3*numel(centres)); P = W; O = W;
for j = 1:numel(centres)
    dz = z-centres(j); s = expm1(1i*dz);
    K = 1i*(s+1)./s; ell = log(abs(s)); t = imag(dz);
    ix = 3*j-2:3*j;
    W(:,ix) = [-2*ell-2i*t.*K, 2i*ell+2*t.*K, 1i*K];
    P(:,ix) = [4*real(K), -4*imag(K), zeros(size(z))];
    O(:,ix) = [-4*imag(K), -4*real(K), zeros(size(z))];
end
end
