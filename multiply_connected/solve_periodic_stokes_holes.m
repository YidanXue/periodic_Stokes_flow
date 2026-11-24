function result = solve_periodic_stokes_holes(geom,dp,disc,options)
% Periodic Stokes flow between graph walls with excluded smooth solids.
% Based on Yidan Xue's poiseuille_smooth.m, commit 16ada0893cdc72b4d20e79a66f4195ca144d83e6.
% Retains VAorthog -> makerows -> least squares -> makefuns architecture.
% Period=2*pi, mu=1, dp=pressure drop per unit length, uv=u+i*v.
% geom: centres, radii, optional walls=[-0.5,0.5], P=2*pi; or h1,h2 as in
%       the original script: h1(x)=x+i*y_top(x), h2(x)=x+i*y_bottom(x).
%       Noncircular solids: geom.obstacles={periodic_ellipse(c,[a,b],angle),...}
%       or structs with boundary(theta), marker, optional derivative(theta).
% disc.holeAAA enables Schwarz-function AAA on each solid (default true for
% geom.obstacles). qHole adds repeated marker poles; zero permits AAA-only.
% disc: nFourier=36, qHole=32, mWall=max(320,16*n), mHole=max(240,12*q).
% options: boundaryVelocity=@(z,component) ..., cutoff=2e-13,
%          completion=true (false is an intentional diagnostic ablation).
% Curved defaults: nFourier=15, at least 600 wall points, wallAAA=true,
% aaaTol=1e-8, aaaMmax=300. Add Chebfun aaa.m to the path for curved walls.
% All evaluators preserve input shape. Evaluate only in fluid/on boundaries.
if nargin<1, geom = struct('centres',pi,'radii',0.2); end
if nargin<2, dp = 1; end
if nargin<3, disc = struct; end
if nargin<4, options = struct; end
geom = defaults(geom,struct('P',2*pi,'walls',[-0.5,0.5],'centres',[],'radii',[]));
geom.centres = geom.centres(:).'; geom.radii = geom.radii(:).';
geom = periodic_obstacles_prepare(geom);
J = numel(geom.centres); geometryInfo = periodic_holes_geometry(geom);
gap = geometryInfo.minimumGap;
validateattributes(dp,{'numeric'},{'real','scalar','finite'});
if geometryInfo.curved, defaultN=15; minWall=600; else, defaultN=36; minWall=320; end
noncircular=isfield(geom,'obstacles');
if noncircular, defaultQ=12; else, defaultQ=32; end
disc = defaults(disc,struct('nFourier',defaultN,'qHole',defaultQ,'holeAAA',noncircular, ...
    'wallAAA',geometryInfo.curved,'aaaTol',1e-8,'aaaMmax',300));
disc = defaults(disc,struct('holeAAATol',disc.aaaTol));
disc.qHole = expand(disc.qHole,J);
disc.holeAAA = expand(disc.holeAAA,J);
assert(all(disc.holeAAA==0 | disc.holeAAA==1),'holeAAA must be logical, scalar or per-obstacle.');
disc = defaults(disc,struct('mWall',max(minWall,16*disc.nFourier),'mHole',max(240,12*disc.qHole)));
disc.mHole = expand(disc.mHole,J);
validateattributes(disc.nFourier,{'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(disc.mWall,{'numeric'},{'scalar','integer','positive','finite'});
if J
    validateattributes(disc.qHole,{'numeric'},{'integer','nonnegative','finite'});
    validateattributes(disc.mHole,{'numeric'},{'integer','positive','finite'});
    assert(all(disc.qHole>0 | disc.holeAAA),'Each obstacle needs centre poles or AAA poles.');
end
options = defaults(options,struct('boundaryVelocity',@(z,c) zeros(size(z)), ...
    'cutoff',2e-13,'completion',true));
validateattributes(options.cutoff,{'numeric'},{'scalar','real','positive','<',1});
[Z,component,weights] = periodic_holes_boundary(geom,disc.mWall,disc.mHole);
Ubc = options.boundaryVelocity(Z,component);
if isscalar(Ubc), Ubc = repmat(Ubc,size(Z)); end
assert(isequal(size(Ubc),size(Z)) && all(isfinite(Ubc)), 'Boundary velocity must be a finite column vector.');
Pol = {}; wallAAA = []; labels = {}; holeAAA=cell(1,J);
if disc.wallAAA || any(disc.holeAAA)
    validateattributes(disc.aaaTol,{'numeric'},{'scalar','real','positive','finite'});
    validateattributes(disc.aaaMmax,{'numeric'},{'scalar','integer','>=',2,'finite'});
end
if any(disc.holeAAA)
    validateattributes(disc.holeAAATol,{'numeric'},{'scalar','real','positive','finite'});
end
if disc.wallAAA
    [Pol,wallAAA] = periodic_wall_aaa(geom,Z(component==1),Z(component==2),disc.aaaTol,disc.aaaMmax);
    for k = 1:numel(wallAAA)
        if ~isempty(wallAAA(k).kept), labels{end+1} = ['wall-',wallAAA(k).source]; end %#ok<AGROW>
    end
end
for j = 1:J
    if disc.holeAAA(j)
        [p,holeAAA{j}] = periodic_obstacle_aaa(periodic_obstacle_at(geom,j), ...
            Z(component==j+2),disc.holeAAATol,disc.aaaMmax);
        if ~isempty(p)
            Pol{end+1}=p; labels{end+1}=sprintf('hole-%d-aaa',j);
        end
        assert(~isempty(p) || disc.qHole(j)>0,'AAA retained no interior poles; add centre terms or refine.');
    end
    if disc.qHole(j)>0
        Pol{end+1} = repmat(geom.centres(j),disc.qHole(j),1);
        labels{end+1} = sprintf('hole-%d-centre',j);
    end
end
[Hes,Rtrain,orthError] = VAorthog(Z,disc.nFourier,Pol);
[R0,R1] = VAeval(Z,Hes,Pol); N = size(R0,2);
[U,V] = makerows(Z,R0,R1);
holeCentres = geom.centres;
if ~options.completion, holeCentres = []; end
[Wh,~,~] = periodic_hole_columns(Z,holeCentres);
A = [U,real(Wh); V,-imag(Wh)];
rhs = [real(Ubc)+dp/2*imag(Z).^2; imag(Ubc)];
map = struct('realF',1:N,'realG',N+(1:N),'a',2*N+1, ...
    'imagF',2*N+1+(1:N),'imagG',3*N+1+(1:N),'dummy',4*N+2, ...
    'hole',4*N+2+(1:3*numel(holeCentres)));
gauges = [N+1,2*N+2,3*N+2,4*N+2];
assert(max(max(abs(A(:,gauges))))==0,'Unexpected nonzero gauge columns.');
active = setdiff(1:size(A,2),gauges);
assert(size(A,1)>numel(active),'Use more boundary samples than active unknowns.');
scale = sqrt(sum(A(:,active).^2,1));
assert(all(isfinite(scale) & scale>0),'Unexpected zero/invalid physical column.');
tic; [Us,Ss,Vs] = svd(A(:,active)./scale,'econ'); sigma = diag(Ss);
keep = sigma>options.cutoff*sigma(1);
c = zeros(size(A,2),1);
c(active) = (Vs(:,keep)*((Us(:,keep)'*rhs)./sigma(keep)))./scale(:);
solveSeconds = toc;
cf = c(map.realF)+1i*c(map.imagF); cg = c(map.realG)+1i*c(map.imagG);
a = c(map.a); ch = reshape(c(map.hole),3,[]);
result.geometry = geom; result.geometryInfo = geometryInfo;
result.dp = dp; result.mu = 1; result.discretization = disc;
result.basis = struct('nBasis',N,'Pol',{Pol},'Hes',{Hes}, ...
    'poleGroupLabels',{labels},'wallAAA',wallAAA,'holeAAA',{holeAAA},'aaaPath',which('aaa'), ...
    'groupOrthogonality',orthError,'recurrenceError',max(abs(Rtrain(:)-R0(:))));
result.boundary = struct('Z',Z,'component',component,'Ubc',Ubc,'quadratureWeights',weights);
result.coefficients = struct('full',c,'cf',cf,'cg',cg,'a',a,'hole',ch,'map',map);
result.solve = struct('rank',sum(keep),'nActive',numel(active),'singularValues',sigma, ...
    'cutoff',options.cutoff,'columnScales',scale,'active',active,'gauges',gauges, ...
    'rowWeighting','none','rows',size(A,1),'seconds',solveSeconds,'minimumGap',gap, ...
    'trainingError',max(abs((A(1:numel(Z),:)*c-rhs(1:numel(Z))) ...
    +1i*(A(numel(Z)+1:end,:)*c-rhs(numel(Z)+1:end)))));
evaluate = @(kind,z) reshape(fh(kind,z(:),cf,cg,a,dp,Hes,Pol,ch,holeCentres),size(z));
result.uv = @(z) evaluate('uv',z); result.p = @(z) evaluate('p',z);
result.omega = @(z) evaluate('omega',z); result.psi = @(z) evaluate('psi',z);
result.inFluid = @(z) periodic_holes_mask(z,geom);
end

function s = defaults(s,d)
names = fieldnames(d);
for j = 1:numel(names)
    if ~isfield(s,names{j}), s.(names{j}) = d.(names{j}); end
end
end

function v = expand(v,J)
if isscalar(v), v = repmat(v,1,J); end
assert(numel(v)==J,'Per-hole parameters must be scalar or have one entry per hole.');
v = v(:).';
end

function [Hes,R,orthError] = VAorthog(Z,n,Pol)
% Original groupwise Vandermonde-Arnoldi, with two MGS passes.
M = numel(Z); Hes = cell(1,2+numel(Pol)); R = []; orthError = zeros(size(Hes));
for group = 1:numel(Hes)
    if group<=2, degree = n; else, degree = numel(Pol{group-2}); end
    Q = ones(M,degree+1); H = zeros(degree+1,degree);
    for k = 1:degree
        if group<=2
            factor = exp(1i*(2*group-3)*Z);
        else
            factor = 1./(exp(1i*Z)-exp(1i*Pol{group-2}(k)));
        end
        q = factor.*Q(:,k);
        initial = norm(q);
        for pass = 1:2
            for j = 1:k
                h = Q(:,j)'*q/M; H(j,k) = H(j,k)+h; q = q-h*Q(:,j);
            end
        end
        H(k+1,k) = norm(q)/sqrt(M);
        assert(H(k+1,k)>100*eps*initial/sqrt(M),'Arnoldi breakdown: reduce degree or change samples.');
        Q(:,k+1) = q/H(k+1,k);
    end
    Hes{group} = H; orthError(group) = norm(Q'*Q/M-eye(degree+1),inf);
    if group==1, R = Q; else, R = [R,Q(:,2:end)]; end %#ok<AGROW>
end
end

function [R0,R1] = VAeval(Z,Hes,Pol)
M = numel(Z); R0 = []; R1 = [];
for group = 1:numel(Hes)
    H = Hes{group}; degree = size(H,2);
    Q = ones(M,degree+1); D = zeros(M,degree+1);
    for k = 1:degree
        if group<=2
            sign = 2*group-3; factor = exp(1i*sign*Z); derivative = 1i*sign*factor;
        else
            ez = exp(1i*Z); factor = 1./(ez-exp(1i*Pol{group-2}(k)));
            derivative = -1i*ez.*factor.^2;
        end
        Q(:,k+1) = (factor.*Q(:,k)-Q(:,1:k)*H(1:k,k))/H(k+1,k);
        D(:,k+1) = (factor.*D(:,k)+derivative.*Q(:,k)-D(:,1:k)*H(1:k,k))/H(k+1,k);
    end
    if group==1, R0 = Q; R1 = D;
    else, R0 = [R0,Q(:,2:end)]; R1 = [R1,D(:,2:end)]; end %#ok<AGROW>
end
end

function [U,V] = makerows(Z,R0,R1)
% Original base packing [Re(cf),Re(cg),a,Im(cf),Im(cg),dummy].
y = imag(Z); zero = zeros(size(Z));
U = [-2*real(R0)+2*y.*imag(R1),real(R1),-4*y, ...
    2*imag(R0)+2*y.*real(R1),-imag(R1),zero];
V = [2*y.*real(R1),-imag(R1),zero,-2*y.*imag(R1),-real(R1),zero];
end

function value = fh(kind,Z,cf,cg,a,dp,Hes,Pol,ch,centres)
[R,D] = VAeval(Z,Hes,Pol); F = R*cf; Fp = D*cf; y = imag(Z);
switch kind
    case 'uv'
        W = -4*a*y-dp/2*y.^2-2*real(F)-2i*y.*Fp+D*cg;
        [Wh,~,~] = periodic_hole_columns(Z,centres);
        value = conj(W+Wh*ch(:));
    case 'p'
        [~,Ph,~] = periodic_hole_columns(Z,centres);
        value = -dp*real(Z)+4*real(Fp)+Ph*ch(:);
    case 'omega'
        [~,~,Oh] = periodic_hole_columns(Z,centres);
        value = 4*a+dp*y-4*imag(Fp)+Oh*ch(:);
    case 'psi'
        value = -2*a*y.^2-dp/6*y.^3-2*y.*real(F)+imag(R*cg);
        for j = 1:numel(centres)
            dz = Z-centres(j); ell = log(abs(expm1(1i*dz)));
            value = value+(ch(3,j)-2*ch(1,j)*imag(dz)).*ell ...
                -2*ch(2,j)*periodic_bloch_wigner(exp(1i*dz));
        end
end
end
