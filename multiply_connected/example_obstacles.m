function r = example_obstacles(shape)
%EXAMPLE_OBSTACLES Solve, check and plot one periodic Stokes flow.
%   r = example_obstacles() uses a disk between flat walls (MATLAB only).
%   'ellipse', 'pair' and 'moving' reproduce Fig. 2(a-c) of the note.
%   These three cases require Chebfun's aaa.m on the MATLAB path.
if nargin == 0
    shape = 'disk';
end

g = struct;
dp = 1;                      % pressure drop per unit length
options.boundaryVelocity = @(z,component) zeros(size(z));
switch shape
    case 'disk'
        g.centres = pi;
        g.radii = 0.2;
        d = struct('nFourier',36,'qHole',32);
    case 'ellipse'
        g.obstacles = {periodic_ellipse(pi,[0.30,0.12],pi/6)};
        d = struct('nFourier',48,'qHole',8,'mHole',600,'holeAAATol',1e-10);
    case {'pair','moving'}
        % Wall functions from poiseuille_smooth.m, PRSA Fig. 3(a).
        g.h1 = @(x) x + 0.5i + 0.2i*sin(2*x);
        g.h2 = @(x) x - 0.5i - 0.2i*sin(2*x);
        d = struct('nFourier',128,'qHole',24,'mHole',800,'holeAAATol',1e-10);
        if strcmp(shape,'pair')
            g.obstacles = {periodic_ellipse(pi-0.5-0.03i,[0.20,0.10],0.35), ...
                           periodic_ellipse(pi+0.5+0.04i,[0.18,0.09],-0.5)};
        else
            particle.marker = pi + 0.25;
            particle.boundary = @(t) particle.marker + 0.21*(1+0.18*cos(3*t)).*exp(1i*t);
            particle.derivative = @(t) 0.21*(-0.54*sin(3*t)+1i*(1+0.18*cos(3*t))).*exp(1i*t);
            g.obstacles = {particle};
            dp = 0;
            options.boundaryVelocity = @(z,component) double(component == 3);
            d.mHole = 1200;
            d.holeAAATol = 1e-11;
        end
    otherwise
        error('Choose disk, ellipse, pair or moving.');
end

r = solve_periodic_stokes_holes(g,dp,d,options);
[z,component] = periodic_holes_boundary(r.geometry,4096,4096,0.371);
boundaryError = max(abs(r.uv(z)-options.boundaryVelocity(z,component)));
fprintf('Maximum velocity error on shifted boundary samples: %.3e\n',boundaryError);
fprintf('Fluid flux at x = 0: %.10f\n',periodic_holes_flux(r,0));
figure('Color','w','Position',[50,50,1200,360]);
plot_periodic_stokes_prsa(r,axes);
end
