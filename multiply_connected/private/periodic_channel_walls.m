function [top,bot] = periodic_channel_walls(g,x)
% Original graph parameterizations h1(x)=x+i*y_top(x), h2(x)=x+i*y_bot(x).
% Evaluate the heights in one period, retaining the caller's horizontal image.
if isfield(g,'h1') && isfield(g,'h2')
    phase = mod(x,2*pi);
    top = x+1i*imag(g.h1(phase)); bot = x+1i*imag(g.h2(phase));
else
    top = x+1i*g.walls(2); bot = x+1i*g.walls(1);
end
end
