function p = plotif(TOPLOT, varargin)
%plot(vargargin) only if TOPLOT is 1
%TOPLOT=1 by default
%EL, 31 July 2016

assert(nargin >= 2);

if isempty(TOPLOT)
    TOPLOT=1;
else
    assert(TOPLOT == 1 || TOPLOT == 0);
end

if TOPLOT==1
    p = plot(varargin{:});
end
    
    
end

