% This code computes the unsteady Couette problem in Fig. 6
%
% Yidan Xue, Jul 2024

warning off
Title = {'(a) $\epsilon=0.1$' '(b) $\epsilon=0.2$' '(c) $\epsilon=0.3$'...
    '(d) $\epsilon=0.5$' '(e) $\epsilon=0.7$' '(f) $\epsilon=0.9$'};
epsilon = [0.1 0.2 0.3 0.5 0.7 0.9];

MS = 'markersize'; LW = 'linewidth'; CO = 'color'; CB = 'colorbar';
FS = 'fontsize'; FW = 'fontweight'; NO = 'normal';
fs = 12;

tiledlayout(3,2,'Padding','tight','TileSpacing','none');
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
polyareac = @(z) polyarea(real(z),imag(z));

for nn = 1:length(epsilon)
    k = 50;
    dt = 2*pi/k;
    tt0 = 0:dt/2:2*pi;
    
    % create function handles for k time steps in a period
    D = 2;
    fre = 1;
    amp = epsilon(nn);
    u_top = 1;
    u_bot = 0;
    dp = 0;   % delta{p}/x
    b = dp/24;

    m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
    h1 = @(x) x + 0.5i*D + 1i*amp*sin(fre*(x));
    h2 = @(x) x - 0.5i*D + 1i*amp*sin(fre*(x));
    top = h1(X); bot = h2(X(end:-1:1));
    top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
    Z = [top;bot];
    Z_bound0 = [top; top(1)+2*pi; bot(end)+2*pi; bot];
    
    uv_list = {};
    
    for jj = 1:k*2
        t = tt0(jj);
        % geometry
        m = 600; X = linspace(0,2*pi,m+1)'; X = X(1:end-1);
        h1 = @(x) x + 0.5i*D + 1i*amp*sin(fre*(x-u_top*t));
        h2 = @(x) x - 0.5i*D + 1i*amp*sin(fre*(x-u_bot*t));
        top = h1(X); bot = h2(X(end:-1:1));
        top_zeta = exp(1i*top); bot_zeta = exp(1i*bot);
        Z = [top;bot];
        Z_bound = [top; top(1)+2*pi; bot(end)+2*pi; bot];
        dp = 0;   % delta{p}/x
        b = dp/24;
    
        % indices
        l1 = 1:m; 
        l2 = m+1:2*m; 
        n = 30; np = 0;  % poly. deg.; # poles per corner
        if amp >= 0.8
            n = 50;
        end
        inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); 
        [r, pol] = aaa(conj(top_zeta),top_zeta);
        pol = -1i*log(pol);
        pol = [pol(real(pol)>=0); pol(real(pol)<0)+2*pi];
        ii = find(~inpolygonc(pol,Z_bound) & real(pol)>=0 & real(pol)<2*pi);
        pol1 = pol(ii);
        [r, pol] = aaa(conj(bot_zeta),bot_zeta);
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
    
        % solution
        c = A\rhs;
        error = A*c-rhs;
        err=max(abs(error));
        [psi,uv,p,omega,f,g] = makefuns(c,Hes,b,Pol);
        uv_list{length(uv_list)+1} = uv;
    end
    
    %% main simulation
    npts = 20; znum = npts*2;
    zpoints1 = pi/2+1i*amp+1i*linspace(-1,1,npts+2).';
    zpoints2 = 3*pi/2-1i*amp+1i*linspace(-1,1,npts+2).';
    zpoints = [zpoints1(2:end-1);zpoints2(2:end-1)];
    zpoints_p = zpoints-1e-8i*sign(imag(zpoints));
    for jj = 1:znum
        if zpoints_p(jj) == zpoints(jj)
            zpoints_p(jj) = zpoints_p(jj)+1e-8;
        end
    end
    
    tt = 0:dt:600*pi; % 600
    for jj = 1:length(tt)-1
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
    
    nexttile
    hold on
    for jj = 1:znum
        if zpoints_mu(jj)>1e-7
            scatter(mod(real(zpoints(jj,:)),2*pi),imag(zpoints(jj,:)),.4,'r')
        else
            scatter(mod(real(zpoints(jj,:)),2*pi),imag(zpoints(jj,:)),.4,'b')
        end
    end
    plot(Z_bound0([1:end 1]),'k',LW,.8)
    title(Title{nn},'interpreter','latex')
    xlim([0 2*pi])
    ylim([-2 2])
    hold off, axis equal off
end

set(gcf,'units','inches','position',[0,0,6,6])
exportgraphics(gcf,'./poincare_section.pdf','Resolution',600)

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