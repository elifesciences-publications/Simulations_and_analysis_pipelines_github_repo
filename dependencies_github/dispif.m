function dispif(TODISP, str)
%disp(str) only if TODISP is 1
%TODISP=1 by default
%EL, 31 July 2016

assert(nargin == 2);

if isempty(TODISP)
    TODISP=1;
else
    assert(TODISP == 1 || TODISP == 0);
end

if TODISP==1
    disp(str);
end
    
    
end

