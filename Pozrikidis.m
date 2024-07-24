% This code computes the steady Couette problem in Fig. 5
%
% Yidan Xue, Jul 2024

MS = 'markersize'; LW = 'linewidth'; CO = 'color'; FS = 'fontsize'; fs = 10;
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
rand(100,100)\rand(100,1);
tic
% parameters
w = 1*pi;
a = 0.2*pi;
dp = 0;   % delta{p}/x
b = -dp/24;
u_top = 1;
u_bot = 0;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 1i*w;
h2 = @(x) x + 1i*a*cos(x);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 25;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(bot_zeta),bot_zeta);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound));
Pol = {pol(ii)};

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
X_p = linspace(-2*pi,4*pi,3*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
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

pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
lev1 = (.1:.1:.9)*pmax;
contour(x,y,pp,lev1,'k',LW,.6)
[xc1,yc1] = contour(x,y,pp,[lev1(2) lev1(2)],'k',LW,.6);
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    y1 = min(imag(Z_p)); y2 = min(xc1(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(200*dy/dmax);
    x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    outside = ~inpolygonc(zz,Z_p);
    pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev2 = (.1:.2:.9)*pmin;
    contour(x,y,pp,lev2,'y',LW,.55)
    [xc2,yc2] = contour(x,y,pp,[1 1]*lev2(2),'y',LW,.55);
    psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
    if sign(fac) == -1     % Moffatt eddies in yellow
        x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
        y1 = min(imag(Z_p)); y2 = min(xc2(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
        dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
        x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
        [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
        outside = ~inpolygonc(zz,Z_p);
        pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
        lev3 = (.1:.2:.9)*pmax;
        contour(x,y,pp,lev3,CO,.99*[1 1 1],LW,.55)
    end
end

plot([cell2mat(Pol)-4*pi cell2mat(Pol)-2*pi cell2mat(Pol) cell2mat(Pol)+2*pi cell2mat(Pol)+4*pi],'.r',MS,5)
ss = sprintf('(a) \\alpha = 0.2\\pi'); title(ss,'FontWeight','Normal')
fontsize(fs,'points')
axis equal off, axis([-2.1*pi 4.1*pi -1.2*pi 1.2*pi]), hold off

tic
% parameters
w = 1*pi;
a = 0.4*pi;
dp = 0;   % delta{p}/x
b = -dp/24;
u_top = 1;
u_bot = 0;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 1i*w;
h2 = @(x) x - 1i*a*cos(x+pi);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 25;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(bot_zeta),bot_zeta);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound));
Pol = {pol(ii)};

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
X_p = linspace(-2*pi,4*pi,3*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
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

pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
lev1 = (.1:.1:.9)*pmax;
contour(x,y,pp,lev1,'k',LW,.6)
[xc1,yc1] = contour(x,y,pp,[lev1(2) lev1(2)],'k',LW,.6);
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    y1 = min(imag(Z_p)); y2 = min(xc1(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(200*dy/dmax);
    x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    outside = ~inpolygonc(zz,Z_p);
    pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev2 = (.1:.2:.9)*pmin;
    contour(x,y,pp,lev2,'y',LW,.55)
    [xc2,yc2] = contour(x,y,pp,[1 1]*lev2(2),'y',LW,.55);
    psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
    if sign(fac) == -1     % Moffatt eddies in yellow
        x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
        y1 = min(imag(Z_p)); y2 = min(xc2(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
        dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
        x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
        [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
        outside = ~inpolygonc(zz,Z_p);
        pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
        lev3 = (.1:.2:.9)*pmax;
        contour(x,y,pp,lev3,CO,.99*[1 1 1],LW,.55)
    end
end

plot([cell2mat(Pol)-4*pi cell2mat(Pol)-2*pi cell2mat(Pol) cell2mat(Pol)+2*pi cell2mat(Pol)+4*pi],'.r',MS,5)
ss = sprintf('(b) \\alpha = 0.4\\pi'); title(ss,'FontWeight','Normal')
fontsize(fs,'points')
axis equal off, axis([-2.1*pi 4.1*pi -1.2*pi 1.2*pi]), hold off

tic
% parameters
w = 1*pi;
a = 0.8*pi;
dp = 0;   % delta{p}/x
b = -dp/24;
u_top = 1;
u_bot = 0;

%setup - geometry
m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
h1 = @(x) x + 1i*w;
h2 = @(x) x - 1i*a*cos(x+pi);
top = h1(X); bot = h2(X(end:-1:1));
top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
Z = [top;bot];
Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];

% indices
l1 = 1:m; 
l2 = m+1:2*m; 
n = 25;
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
[r, pol] = aaa(conj(bot_zeta),bot_zeta);
pol = -1i*log(pol);
pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
ii = find(~inpolygonc(pol,Z_bound));
Pol = {pol(ii)};

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
X_p = linspace(-2*pi,4*pi,3*m+1)';
top_p = h1(X_p(end:-1:1)); bot_p = h2(X_p); Z_p = [top_p;bot_p];
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

pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
lev1 = (.2:.2:.8)*pmax;
contour(x,y,pp,lev1,'k',LW,.6)
[xc1,yc1] = contour(x,y,pp,[lev1(2) lev1(2)],'k',LW,.6);
psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
if sign(fac) == -1     % Moffatt eddies in yellow
    x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    y1 = min(imag(Z_p)); y2 = min(xc1(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    dmax = max(dx,dy); nx = ceil(600*dx/dmax); ny = ceil(200*dy/dmax);
    x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    outside = ~inpolygonc(zz,Z_p);
    pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    lev2 = (.1:.2:.9)*pmin;
    contour(x,y,pp,lev2,'y',LW,.55)
    [xc2,yc2] = contour(x,y,pp,[1 1]*lev2(2),'y',LW,.55);
    psiratio = pmin/pmax; fac = max(psiratio,1/psiratio);
    if sign(fac) == -1     % Moffatt eddies in yellow
        x1 = min(real(Z_p)); x2 = max(real(Z_p)); xm = mean([x1 x2]); dx = diff([x1 x2]);
        y1 = min(imag(Z_p)); y2 = min(xc2(2,:)); ym = mean([y1 y2]); dy = diff([y1 y2]);
        dmax = max(dx,dy); nx = ceil(400*dx/dmax); ny = ceil(400*dy/dmax);
        x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
        [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
        outside = ~inpolygonc(zz,Z_p);
        pp = psi(zz)-psi(-pi/2); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
        lev3 = (.1:.2:.9)*pmax;
        contour(x,y,pp,lev3,CO,.99*[1 1 1],LW,.55)
    end
end

plot([cell2mat(Pol)-4*pi cell2mat(Pol)-2*pi cell2mat(Pol) cell2mat(Pol)+2*pi cell2mat(Pol)+4*pi],'.r',MS,5)
ss = sprintf('(c) \\alpha = 0.8\\pi'); title(ss,'FontWeight','Normal')
fontsize(fs,'points')
axis equal off, axis([-2.1*pi 4.1*pi -1.2*pi 1.2*pi]), hold off

exportgraphics(gcf,'Pozrikidis.pdf','Resolution',600)

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
% cZ = spdiags(conj(Z),0,M,M);                    % conj(Z)
ZZ = spdiags(Z,0,M,M);                    
% PSI = [-2*imag(ZZ)*real(R0) imag(R0) -2*imag(Z).^2 2*imag(ZZ)*imag(R0) real(R0) zeros(M,1)];
U = [-2*real(R0)+2*imag(ZZ)*imag(R1) real(R1) -4*imag(Z) 2*imag(R0)+2*imag(ZZ)*real(R1) -imag(R1) zeros(M,1)];
V = [2*imag(ZZ)*real(R1) -imag(R1) zeros(M,1) -2*imag(ZZ)*imag(R1) -real(R1) zeros(M,1)];
% P = [4*real(R1) zero zeros(M,1) -4*imag(R1) zero zeros(M,1)];
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