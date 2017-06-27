%2016-08-23, EL: return current date and time as a string, useful to make
%unique file names
%allow to pass an argument to specify format of date;
%by default, without an argument, spit out date as 'yyyy-mm-dd_HH.MM.SS'
function out = getDate(varargin)
    if nargin < 1 
        out = datestr(now,'yyyy-mm-dd_HH.MM.SS');
        %disp(['nargin ' num2str(nargin)]);
    else
        assert(ischar(varargin{1}))
        out = datestr(now,varargin{1});
    end
end
%getDate = @(x) datestr(now,'yyyy-mm-dd_HH.MM.SS');