% This code computes the unsteady Couette problem in Fig. 9c
%
% Yidan Xue, Jul 2024

warning off
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
polyareac = @(z) polyarea(real(z),imag(z));

MS = 'markersize'; LW = 'linewidth'; CO = 'color'; CB = 'colorbar';
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal';
fs = 12;

k = 50; % 50
dt = 2*pi/k;
tt0 = 0:dt/2:2*pi;

% create function handles for k time steps in a period
D = 1;
amp = 1;
u_top = 1;
u_bot = 0;
dp = 0;   % delta{p}/x
b = dp/24;

w1 = pi/2+D/2*1i; w2 = pi/2+(D/2+amp)*1i; w3 = 3*pi/2+(D/2+amp)*1i; w4 = 3*pi/2+D/2*1i;
w = [w1;w2;w3;w4];
m = 600; s = tanh(linspace(-14,14,m));
s1 = tanh(linspace(-14,0,m/2)); s2 = tanh(linspace(0,14,m/2));
top = [(w4-2*pi+w1)/2+(w1-w4+2*pi)/2*s2 (w2+w1)/2+(w2-w1)/2*s (w3+w2)/2+(w3-w2)/2*s...
    (w4+w3)/2+(w4-w3)/2*s (w1+2*pi+w4)/2+(w1+2*pi-w4)/2*s1].';
bot = flip(conj(top));
Z = [top; bot]; 
% w = [w;flip(w-D*1i)];
Z_bound0 = [top; top(1)+2*pi; bot(end)+2*pi; bot];
uv_list = {};

for jj = 1:k*2
    t = tt0(jj);
    %setup - geometry
    top_t = top+u_top*t;
    top_t = [top_t(real(top_t)>2*pi)-2*pi; top_t(real(top_t)<=2*pi)];
    Z = [top_t; bot]; 
    w_t = w+u_top*t;
    w_b = flip(conj(w));
    w_t = [w_t(real(w_t)>2*pi)-2*pi; w_t(real(w_t)<=2*pi)];
    Z_bound = [top_t-2*pi; top_t; top_t+2*pi; bot+2*pi; bot; bot-2*pi];
    zeta_t = exp(1i*top_t); zeta_b = exp(1i*bot);
    zeta_wt = exp(1i*w_t); zeta_wb = exp(1i*w_b); 
    
    % indices
    l1 = 1:4*m; 
    l2 = 4*m+1:8*m; 
    n = 30; % poly. deg.;
    np = 24;                                    % poles per corner
    l = 1; dk = l*cluster(np,4);
    t1 = 3*pi/4; t2 = 3*pi/4; t3 = pi/4; t4 = pi/4;
    pol1 = w(1)+exp(1i*t1)*dk; pol2 = w(2)+exp(1i*t2)*dk;
    pol3 = w(3)+exp(1i*t3)*dk; pol4 = w(4)+exp(1i*t4)*dk;
    Pol = {pol1+u_top*t,pol2+u_top*t,pol3+u_top*t,pol4+u_top*t,...
        conj(pol1),conj(pol2),conj(pol3),conj(pol4)};
    % Pol = {};
    % for nk = 1:length(zeta_wt)
    %     ii = find(abs(zeta_t-zeta_wt(nk)) == min(abs(zeta_t-zeta_wt.'),[],2));
    %     [r,pol] = aaa(conj(zeta_t(ii)),zeta_t(ii),'tol',1e-8);
    %     pol = -1i*log(pol);
    %     pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    %     jj = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
    %     Pol{length(Pol)+1} = pol(jj).';
    % end
    % for nk = 1:length(zeta_wb)
    %     ii = find(abs(zeta_b-zeta_wb(nk)) == min(abs(zeta_b-zeta_wb.'),[],2));
    %     [r,pol] = aaa(conj(zeta_b(ii)),zeta_b(ii),'tol',1e-8);
    %     pol = -1i*log(pol);
    %     pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    %     jj = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
    %     Pol{length(Pol)+1} = pol(jj).';
    % end
    % [r, pol] = aaa(conj(top_zeta),top_zeta, 'mmax', 200);
    % pol = -1i*log(pol);
    % pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    % ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
    % pol1 = pol(ii);
    % [r, pol] = aaa(conj(bot_zeta),bot_zeta, 'mmax', 200);
    % pol = -1i*log(pol);
    % pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
    % ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
    % pol2 = pol(ii);
    % Pol = {pol1;pol2};
    
    Hes = VAorthog(Z,n,Pol);   % Arnoldi Hessenberg matrices
    
    % boundary conditions
    
    [A1,rhs1,A2,rhs2,U,V] = makerows(Z,n,Hes,Pol);   
    A1(l1,:) =   U(l1,:); rhs1(l1) = u_top+12*b*imag(Z(l1)).^2;  
    A2(l1,:) =   V(l1,:); rhs2(l1) = 0;
    A1(l2,:) =   U(l2,:); rhs1(l2) = u_bot+12*b*imag(Z(l2)).^2;  
    A2(l2,:) =   V(l2,:); rhs2(l2) = 0;
    A = [A1; A2]; rhs = [rhs1; rhs2];
    
    % solution and plot
    [A,rhs] = rowweighting(A,rhs,Z,[w_t;w_b]);
    c = A\rhs;
    error = A*c-rhs;
    % semilogy(abs(error),'.')
    err=max(abs(error))
    [psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);
    % x1 = min(real(Z_bound)); x2 = max(real(Z_bound)); xm = mean([x1 x2]); dx = diff([x1 x2]);
    % y1 = min(imag(Z_bound)); y2 = max(imag(Z_bound)); ym = mean([y1 y2]); dy = diff([y1 y2]);
    % dmax = max(dx,dy); nx = ceil(300*dx/dmax); ny = ceil(300*dy/dmax);
    % x = linspace(x1,x2,nx); y = linspace(y1,y2,ny);
    % [xx,yy] = meshgrid(x,y); zz = xx + 1i*yy;
    % inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
    % outside = ~inpolygonc(zz,Z_bound);
    % uu = abs(uv(zz)); uu(outside) = NaN; umax = max(max(uu));
    % pcolor(x,y,uu), hold on, colormap(gca,parula)
    % shading interp, c=colorbar, caxis([0 umax])
    % c.Label.String = 'Velocity magnitude';
    % plot(Z_bound([1:end 1]),'k',LW,.8), hold on
    % pp = psi(zz); pp(outside) = NaN; pmin = min(min(pp)); pmax = max(max(pp));
    % lev = pmin+(.1:.1:.9)*(pmax-pmin);
    % contour(x,y,pp,lev,'k',LW,.6)
    % plot(cell2mat(Pol),'.r',MS,8)
    % fontsize(fs,'points')
    % axis equal off, axis([0 2*pi -2 2]), hold off
    uv_list{length(uv_list)+1} = uv;
    % pause
end
    
%% main simulation
npts = 20; znum = npts;
% zpoints1 = (1-amp)*1i-1i*D*linspace(0,1,npts+2).';
% zpoints2 = pi/2+1i-1i*D*linspace(0,1,npts+2).';
% zpoints3 = pi+(1-amp)*1i-1i*D*linspace(0,1,npts+2).';
% zpoints4 = 3*pi/2+1i-1i*D*linspace(0,1,npts+2).';
% zpoints2 = 3*pi/2-1i*amp+1i*linspace(-1,1,npts+2).';
% zpoints3 = pi+1i*linspace(-1,1,npts+2).';
% zpoints = [zpoints1(2:end-1); zpoints2(2:end-1); zpoints3(2:end-1); zpoints4(2:end-1)];
zpoints1 = pi+1i*(D/2+amp)*linspace(-1,1,npts+2).';
zpoints = zpoints1(2:end-1);
% zpoints = [zpoints; 5*pi/8-0.1i; 5*pi/8; 5*pi/8+0.1i; 11*pi/8-0.1i; 11*pi/8; 11*pi/8+0.1i];
zpoints_p = zpoints-1e-8i*sign(imag(zpoints));
for jj = 1:znum
    if zpoints_p(jj) == zpoints(jj)
        zpoints_p(jj) = zpoints_p(jj)+1e-8;
    end
end

%%
tt = 0:dt:600*pi; % 600
for jj = 1:length(tt)-1
    if mod(jj,k*50) == 0
        jj/(length(tt)-1)
    end
    num = mod(jj-1,k);
    uv1 = uv_list{mod(2*num,2*k)+1};
    uv2 = uv_list{mod(2*num+1,2*k)+1};
    uv3 = uv_list{mod(2*num+2,2*k)+1};
    % 4th-order Runge-Kutta
    pts = zpoints(:,end);
    k1 = uv1(pts);
    k2 = uv2(pts+dt*k1/2);
    k3 = uv2(pts+dt*k2/2);
    k4 = uv3(pts+dt*k3);
    pts_new = pts+dt/6*(k1+2*k2+2*k3+k4);
    zpoints = [zpoints pts_new];
    pts = zpoints_p(:,end);
    k1 = uv1(pts);
    k2 = uv2(pts+dt*k1/2);
    k3 = uv2(pts+dt*k2/2);
    k4 = uv3(pts+dt*k3);
    pts_new = pts+dt/6*(k1+2*k2+2*k3+k4);
    zpoints_p = [zpoints_p pts_new];
end

zpoints = zpoints(:,1:k:end);
zpoints_p = zpoints_p(:,1:k:end);
zpoints_mu = zeros(znum,1);
for jj = 1:znum
    zpoints_mu(jj) = sum(abs(zpoints(jj,2:101)-zpoints_p(jj,2:101)).^2);
end
%%
hold on
for jj = 1:znum
    if zpoints_mu(jj)>1e-7
    % if abs(zpoints(jj,100)-zpoints_p(jj,100))>1e-4
    % if max(abs(zpoints(jj,1:100)-zpoints_p(jj,1:100)))>1e-4
        scatter(mod(real(zpoints(jj,:)),2*pi),imag(zpoints(jj,:)),.4,'r')
    else
        scatter(mod(real(zpoints(jj,:)),2*pi),imag(zpoints(jj,:)),.4,'b')
    end
end
plot(Z_bound0([1:end 1]),'k',LW,.8)
xlim([0 2*pi])
ylim([-1.5 1.5])
hold off, axis equal off

set(gcf,'units','inches','position',[0,0,6,6])
exportgraphics(gcf,'./poincare_case3.pdf','Resolution',600)

function d = cluster(n,sigma)                     
    nc = ceil(n); d = exp(sigma*(sqrt(nc:-1:1)-sqrt(nc)));
end

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

function [A,rhs] = rowweighting(A,rhs,Z,w)
dZw = min(abs(Z-w.'),[],2);
wt = [dZw; dZw];
M2 = 2*length(Z); W = spdiags(wt,0,M2,M2);
A = W*A; rhs = W*rhs;
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