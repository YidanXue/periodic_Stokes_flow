% This code computes the Poiseuille problem in Fig. 4
%
% Yidan Xue, Jul 2024

% case a
MS = 'markersize'; LW = 'linewidth'; CO = 'color'; FS = 'fontsize'; fs = 12;
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
tic
% parameters
D = 1;
u_top = 0;
u_bot = 0;
dp = 1;   % delta{p}/x
b = dp/24;

%setup - geometry
w1 = pi/3+0.5i*D; w2 = 2*pi/3+1i*D; w3 = 4*pi/3+1i*D; w4 = 5*pi/3+0.5i*D;
w = [w1;w2;w3;w4]; w = [w;conj(w)];
m = 300; s = tanh(linspace(-7,7,m));
s1 = tanh(linspace(-7,0,m/2)); s2 = tanh(linspace(0,7,m/2));
top = [(w4-2*pi+w1)/2+(w1-w4+2*pi)/2*s2 (w2+w1)/2+(w2-w1)/2*s (w3+w2)/2+(w3-w2)/2*s...
    (w4+w3)/2+(w4-w3)/2*s (w1+2*pi+w4)/2+(w1+2*pi-w4)/2*s1].';
bot = flip(conj(top));
Z = [top; bot]; 
w_zeta = exp(1i*w(real(w)<=pi));
zeta = exp(1i*Z(real(Z)<=pi)); 

% indices
l1 = 1:4*m; 
l2 = 4*m+1:8*m; 
n = 15; % poly. deg.;
Pol = {};
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
for k = 1:length(w_zeta)
    ii = find(abs(zeta-w_zeta(k)) == min(abs(zeta-w_zeta.'),[],2));
    [r,pol] = aaa(conj(zeta(ii)),zeta(ii),'tol',1e-8);
    pol = -1i*log(pol);
    pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    jj = find(~inpolygonc(pol,Z) & real(pol)>=0 & real(pol)<2*pi);
    Pol{length(Pol)+1} = pol(jj).';
    Pol{length(Pol)+1} = 2*pi-conj(pol(jj)).';
end

Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices

% boundary conditions

[A1,rhs1,A2,rhs2,U,V] = makerows(Z,n,Hes,Pol);   
A1(l1,:) =   U(l1,:); rhs1(l1) = u_top+12*b*imag(Z(l1)).^2;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
A1(l2,:) =   U(l2,:); rhs1(l2) = u_bot+12*b*imag(Z(l2)).^2;  
A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
c = A\rhs;
tsolve=toc
error = A*c-rhs;
err=max(abs(error));
[psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);

nexttile
Z_p = [1i*D-pi; w3-2*pi; w4-2*pi; w1; w2; w3; w4; w1+2*pi; w2+2*pi; 1i*D+3*pi];
Z_p = [Z_p;flip(conj(Z_p))];
x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z_p)); y2 = max(imag(Z_p)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(600*dy/dmax);
x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
outside = ~inpolygonc(zz,Z_p);
uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
pcolor(x,y,uu), hold on, colormap(gca,parula)
shading interp, c=colorbar, caxis([0 umax])
c.Label.String = 'Velocity magnitude';
plot(Z_p([1:end 1]),'k',LW,.8), hold on
pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
lev = pmin+(.1:.1:.9)*(pmax-pmin);
contour(x,y,pp,lev,'k',LW,.6)
plot([cell2mat(Pol)-2*pi cell2mat(Pol) cell2mat(Pol)+2*pi],'.r',MS,8)
plot([0,0],[-3,3],'k--',LW,1.2)
plot([2*pi,2*pi],[-3,3],'k--',LW,1.2)
text(-pi+0.3,1.7,'(a)','horizontalalignment', 'center','fontsize',12)
fontsize(fs,'points')
axis equal off, axis([-pi 3*pi -2 2]), hold off

%%
tic
% parameters
D = 1;
u_top = 0;
u_bot = 0;
dp = 1;   % delta{p}/x
b = dp/24;

%setup - geometry
w1 = pi/3; w2 = pi/3+1i*D; w3 = 5*pi/3+1i*D; w4 = 5*pi/3;
w5 = 4*pi/3-1i*D; w6 = 4*pi/3; w7 = 2*pi/3; w8 = 2*pi/3-1i*D;
w = [w1;w2;w3;w4;w5;w6;w7;w8];
m = 400; s = tanh(linspace(-7,7,m));
s1 = tanh(linspace(-7,0,m/2)); s2 = tanh(linspace(0,7,m/2));
top = [(w4-2*pi+w1)/2+(w1-w4+2*pi)/2*s2 (w2+w1)/2+(w2-w1)/2*s (w3+w2)/2+(w3-w2)/2*s...
    (w4+w3)/2+(w4-w3)/2*s (w1+2*pi+w4)/2+(w1+2*pi-w4)/2*s1].';
bot = [(w5+w8+2*pi)/2+(w5-w8-2*pi)/2*s2 (w6+w5)/2+(w6-w5)/2*s (w7+w6)/2+(w7-w6)/2*s...
    (w8+w7)/2+(w8-w7)/2*s (w5-2*pi+w8)/2+(w5-2*pi-w8)/2*s1].';
Z = [top; bot];
w_zeta = exp(1i*w(real(w)<=pi));
zeta = exp(1i*Z(real(Z)<=pi)); 

% indices
l1 = 1:4*m; 
l2 = 4*m+1:8*m; 
n = 15; % poly. deg.;
Pol = {};
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
for k = 1:length(w_zeta)
    ii = find(abs(zeta-w_zeta(k)) == min(abs(zeta-w_zeta.'),[],2));
    [r,pol] = aaa(conj(zeta(ii)),zeta(ii),'tol',1e-8);
    pol = -1i*log(pol);
    pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    jj = find(~inpolygonc(pol,Z) & real(pol)>=0 & real(pol)<2*pi);
    Pol{length(Pol)+1} = pol(jj).';
    Pol{length(Pol)+1} = 2*pi-conj(pol(jj)).';
end

Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices

% boundary conditions

[A1,rhs1,A2,rhs2,U,V] = makerows(Z,n,Hes,Pol);   
A1(l1,:) =   U(l1,:); rhs1(l1) = u_top+12*b*imag(Z(l1)).^2;  
A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
A1(l2,:) =   U(l2,:); rhs1(l2) = u_bot+12*b*imag(Z(l2)).^2;  
A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
A = [A1; A2]; rhs = [rhs1; rhs2];

% solution and plot
c = A\rhs;
tsolve=toc
error = A*c-rhs;
err=max(abs(error));
[psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);

nexttile
Z_p = [-pi+1i*D; w3-2*pi; w4-2*pi; w1; w2; w3; w4; w1+2*pi; w2+2*pi; 3*pi+1i*D;...
    3*pi; w7+2*pi; w8+2*pi; w5; w6; w7; w8; w5-2*pi; w6-2*pi; -pi];
x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z_p)); y2 = max(imag(Z_p)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(600*dy/dmax);
x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
outside = ~inpolygonc(zz,Z_p);
uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
pcolor(x,y,uu), hold on, colormap(gca,parula)
shading interp, c=colorbar, caxis([0 umax])
c.Label.String = 'Velocity magnitude';
plot(Z_p([1:end 1]),'k',LW,.8), hold on
pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
lev = pmin+(.1:.1:.9)*(pmax-pmin);
contour(x,y,pp,lev,'k',LW,.6)
plot([cell2mat(Pol)-2*pi cell2mat(Pol) cell2mat(Pol)+2*pi],'.r',MS,8)
plot([0,0],[-3,3],'k--',LW,1.2)
plot([2*pi,2*pi],[-3,3],'k--',LW,1.2)

% top eddies
x = linspace(real(w2)-0.5,real(w2)+0.5,300); y = linspace(imag(w2)-0.5,imag(w2)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi+1i*D); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmax;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w3)-0.5,real(w3)+0.5,300); y = linspace(imag(w3)-0.5,imag(w3)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi+1i*D); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmax;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w2)+2*pi-0.5,real(w2)+2*pi+0.5,300); y = linspace(imag(w2)-0.5,imag(w2)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi+1i*D); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmax;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w3)-2*pi-0.5,real(w3)-2*pi+0.5,300); y = linspace(imag(w3)-0.5,imag(w3)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi+1i*D); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmax;
    contour(x,y,pp,lev,'y',LW,.55)
end

% bottom eddies
x = linspace(real(w5)-0.5,real(w5)+0.5,300); y = linspace(imag(w5)-0.5,imag(w5)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmin;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w8)-0.5,real(w8)+0.5,300); y = linspace(imag(w8)-0.5,imag(w8)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmin;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w5)-2*pi-0.5,real(w5)-2*pi+0.5,300); y = linspace(imag(w5)-0.5,imag(w5)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmin;
    contour(x,y,pp,lev,'y',LW,.55)
end

x = linspace(real(w8)+2*pi-0.5,real(w8)+2*pi+0.5,300); y = linspace(imag(w8)-0.5,imag(w8)+0.5,300);
[xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
outside = ~inpolygonc(zz,Z_p);
pp = psi(zz)-psi(pi); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    lev = (.1:.4:.9)*pmin;
    contour(x,y,pp,lev,'y',LW,.55)
end

text(-pi+0.3,1.7,'(b)','horizontalalignment', 'center','fontsize',12)
fontsize(fs,'points')
axis equal off, axis([-pi 3*pi -2 2]), hold off
exportgraphics(gcf,'poiseuille_sharp.pdf','Resolution',600)

function [Hes,R] = VAorthog(Z,n,varargin)  % Vand.+Arnoldi orthogonalization
M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end
% First orthogonalize the polynomial part
Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
   q = exp(-1i*Z).*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{1} = H; R = Q;

Q = ones(M,1); H = zeros(n+1,n);
for k = 1:n
   q = exp(1i*Z).*Q(:,k);
   for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
   H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
end
Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];

% Next orthogonalize the pole parts, if any
while ~isempty(Pol)
   pol = Pol{1}; Pol(1) = [];
   np = length(pol); H = zeros(np,np-1); Q = ones(M,1);
   for k = 1:np
      q = Q(:,k)./(exp(1i*Z)-exp(1i*pol(k)));
      for j = 1:k, H(j,k) = Q(:,j)'*q/M; q = q - H(j,k)*Q(:,j); end
      H(k+1,k) = norm(q)/sqrt(M); Q(:,k+1) = q/H(k+1,k);
    end
   Hes{length(Hes)+1} = H; R = [R Q(:,2:end)];
end
end

function [R0,R1] = VAeval(Z,Hes,varargin)  % Vand.+Arnoldi basis construction

M = length(Z); Pol = []; if nargin == 3, Pol = varargin{1}; end

     H = Hes{1}; Hes(1) = []; n = size(H,2);
     Q = ones(M,1); D = zeros(M,1);
     Zpki = exp(-1i*Z); Zpkid = -1i*exp(-1i*Z);
    for k = 1:n
        hkk = H(k+1,k);
        Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
        D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
    end
     R0 = Q; R1 = D;

     H = Hes{1}; Hes(1) = []; n = size(H,2);
     Q = ones(M,1); D = zeros(M,1);
     Zpki = exp(1i*Z); Zpkid = 1i*exp(1i*Z);
    for k = 1:n
        hkk = H(k+1,k);
        Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
        D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
    end
     R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
    
     % Next construct the pole parts of the basis, if any
     while ~isempty(Pol)
        pol = Pol{1}; Pol(1) = [];
        H = Hes{1}; Hes(1) = []; np = length(pol); Q = ones(M,1); D = zeros(M,1);
        for k = 1:np
            Zpki = 1./(exp(1i*Z)-exp(1i*pol(k))); Zpkid = -exp(1i*Z)*1i./(exp(1i*Z)-exp(1i*pol(k))).^2;
            hkk = H(k+1,k);
            Q(:,k+1) = ( Q(:,k).*Zpki - Q(:,1:k)*H(1:k,k))/hkk;
            D(:,k+1) = ( D(:,k).*Zpki - D(:,1:k)*H(1:k,k) + Q(:,k).*Zpkid )/hkk;
        end
        R0 = [R0 Q(:,2:end)]; R1 = [R1 D(:,2:end)];
     end
end

function [A1,rhs1,A2,rhs2,U,V] = makerows(Z,n,Hes,varargin)
Pol = []; if nargin == 4, Pol = varargin{1}; end
[R0,R1] = VAeval(Z,Hes,Pol);
M = length(Z); N = 4*size(R0,2)+2; zero = 0*R0;
ZZ = spdiags(Z,0,M,M);                    
U = [-2*real(R0)+2*imag(ZZ)*imag(R1) real(R1) -4*imag(Z) 2*imag(R0)+2*imag(ZZ)*real(R1) -imag(R1) zeros(M,1)];
V = [2*imag(ZZ)*real(R1) -imag(R1) zeros(M,1) -2*imag(ZZ)*imag(R1) -real(R1) zeros(M,1)];
A1 = zeros(M,N); rhs1 = zeros(M,1);
A2 = zeros(M,N); rhs2 = zeros(M,1);
end

% function [A,rhs] = rowweighting(A,rhs,Z,w)
% dZw = min(abs(Z-w.'),[],2);
% wt = [dZw; dZw];
% M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
% A = W*A; rhs = W*rhs;
% end

function [psi,uv,p,omega,f,g] = makefuns(c,Hes,b,varargin)  % make function handles
Pol = []; if nargin == 4, Pol = varargin{1}; end
cc = c(1:end/2) + 1i*c(end/2+1:end);
reshaper = @(str) @(z) reshape(fh(str,z(:),cc,Hes,b,Pol),size(z));
  psi = reshaper('psi');    uv = reshaper('uv');    p = reshaper('p');
omega = reshaper('omega');   f = reshaper('f');   g = reshaper('g');
end

function fh = fh(i,Z,cc,Hes,b,Pol)
[R0,R1] = VAeval(Z,Hes,Pol);
N = size(R0,2);
cf = cc(1:N); cg = cc(N+(1:N)); a = real(cc(2*N+1));
ff = R0*cf; gg = R0*cg;
switch i
   case   'f'  , fh = -1i*a*Z-3*b*Z.^2+ff;
   case   'g'  , fh = 1i*a*Z.^2+b*Z.^3-Z.*ff+gg;
   case  'psi' , fh = -2*a*imag(Z).^2-4*b*imag(Z).^3-2*imag(Z).*real(ff)+imag(gg);
   case   'uv' , fh = conj(-4*a*imag(Z)-12*b*imag(Z).^2-2*real(ff)-2i*imag(Z).*(R1*cf)+R1*cg);
   case   'p'  , fh = -24*b*real(Z)+4*real(R1*cf);                      
   case 'omega', fh = 4*a+24*b*imag(Z)-4*real(R1*cf);                       
end
end