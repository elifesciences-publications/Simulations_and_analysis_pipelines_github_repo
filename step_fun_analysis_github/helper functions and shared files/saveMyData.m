function saveMyData(VAR, FILENAME, OUTDIR)
%Save variable VAR as an ascii file FILENAME.txt in OUTDIR.
%   Filename gets timestamped. FILENAME, OUTDIR optional, but must be
%   strings if passed. VARNAME must be a string.

if ~isempty(FILENAME)
    assert(ischar(FILENAME));
end
if ~isempty(OUTDIR)
    assert(ischar(OUTDIR));
else
    OUTDIR='.';
end

%get name of variable
varname = @(x) inputname(1);
VARNAME = varname(VAR);

%save
TIMEFORMAT = 'yyyy-mm-dd_HH.MM.SS';
TIMESTAMP = datestr(now,TIMEFORMAT);
save([OUTDIR '/' TIMESTAMP '_' FILENAME '.mat'], VARNAME, '-mat');

end

