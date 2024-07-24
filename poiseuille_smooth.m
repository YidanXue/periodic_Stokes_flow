% This code computes the Poiseuille problem in Fig. 3
%
% Yidan Xue, Jul 2024

% case a
MS = 'markersize'; LW = 'linewidth'; CO = 'color'; FS = 'fontsize'; fs = 12;
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
tic
% parameters
D = 1;
fre1 = 2;
fre2 = 2;
amp1 = 0.2;
amp2 = 0.2;
u_top = 0;
u_bot = 0;
dp = 1;   % delta{p}/x
b = dp/24;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 0.5i*D + 1i*amp1*sin(fre1*x);
h2 = @(x) x - 0.5i*D - 1i*amp2*sin(fre2*x);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];
 
% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 15; np = 0;  % poly. deg.; # poles per corner
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(top_zeta),top_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol1 = pol(ii);
[r, pol] = aaa(conj(bot_zeta),bot_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol2 = pol(ii);
Pol = {pol1;pol2};

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
tsolve=toc;
error = A*c-rhs;
err=max(abs(error));
[psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);

nexttile
X_p = linspace(-pi,3*pi,2*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z_p)); y2 = max(imag(Z_p)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(800*dx/dmax); ny = ceil(800*dy/dmax);
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
text(-pi+0.3,1.3,'(a)','horizontalalignment', 'center','fontsize',12)
fontsize(fs,'points')
axis equal off, axis([-pi 3*pi -1.6 1.6]), hold off

%% case b
tic
% parameters
D = 1;
fre1 = 1;
fre2 = 2;
amp1 = 0.2;
amp2 = 0.3;
u_top = 0;
u_bot = 0;
dp = 1;   % delta{p}/x
b = dp/24;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 0.5i*D + 1i*amp1*sin(fre1*(x+pi/4));
h2 = @(x) x - 0.5i*D - 1i*amp2*sin(fre2*x);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];
 
% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 15; np = 0;  % poly. deg.; # poles per corner
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(top_zeta),top_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol1 = pol(ii);
[r, pol] = aaa(conj(bot_zeta),bot_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol2 = pol(ii);
Pol = {pol1;pol2};

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
tsolve=toc;
error = A*c-rhs;
err=max(abs(error));
[psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);

nexttile
X_p = linspace(-pi,3*pi,2*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z_p)); y2 = max(imag(Z_p)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(800*dx/dmax); ny = ceil(800*dy/dmax);
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
text(-pi+0.3,1.3,'(b)','horizontalalignment', 'center','fontsize',12)
fontsize(fs,'points')
axis equal off, axis([-pi 3*pi -1.6 1.6]), hold off

%% case c
tic
% parameters
D = 1;
fre1 = 1;
fre2 = 1;
amp1 = .5;
amp2 = .5;
u_top = 0;
u_bot = 0;
dp = 1;   % delta{p}/x
b = dp/24;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 0.5i*D + 1i*amp1*(tanh(cos(1+2*sin(x)).^2)-0.5);
h2 = @(x) x - 0.5i*D - 1i*amp2*(tanh(cos(1+2*sin(x)).^2)-0.5);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];
 
% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 15; np = 0;  % poly. deg.; # poles per corner
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(top_zeta),top_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol1 = pol(ii);
[r, pol] = aaa(conj(bot_zeta),bot_zeta,'tol',1e-8,'mmax',300);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
pol2 = pol(ii);
Pol = {pol1;pol2};

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
tsolve=toc;
error = A*c-rhs;
err=max(abs(error));
[psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);

nexttile
X_p = linspace(-pi,3*pi,2*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
y1 = min(imag(Z_p)); y2 = max(imag(Z_p)); ym = mean([y1 y2]); dy = diff([y1 y2]);
dmax = max(dx,dy); nx = ceil(800*dx/dmax); ny = ceil(800*dy/dmax);
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
text(-pi+0.3,1.3,'(c)','horizontalalignment', 'center','fontsize',12)
fontsize(fs,'points')
axis equal off, axis([-pi 3*pi -1.6 1.6]), hold off

exportgraphics(gcf,'poiseuille_smooth.pdf','Resolution',600)

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