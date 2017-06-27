function [dsFinal] = joinSheets(INFILE, SHEETS)
%%2016-02-13, EL
%%Script to join different sheets of an Excel file from the Tecan plate
%%reader. Use timestamps on each sheet to compute the relative time of
%%measurements. Useful when measurement gets interrupted and restarted.
%%Pass in INFILE as string and SHEETS as a cell array of strings that label
%%Excel sheets. 

dsFinal = dataset;

for snum=1:numel(SHEETS)

%get time and date
[dnum,dtxt,draw] = xlsread(INFILE,SHEETS{snum},'B5');
[tnum,ttxt,traw] = xlsread(INFILE,SHEETS{snum},'B6');

%fix them
[fvec{snum}, fstr{snum}, fserial{snum}] = fixExcelTime(datenum(dnum), traw);

%get measurement interval time
[mnum, mtxt, mraw] = xlsread(INFILE, SHEETS{snum}, 'E27');

%column A contains names of parameters; column E contains values
[Anum, Atxt, Araw] = xlsread(INFILE, SHEETS{snum}, 'A:A');
[Enum, Etxt, Eraw] = xlsread(INFILE, SHEETS{snum}, 'E:E');

%find 'Kinetic Cyles'
ixKin=find(cellfun(@(s) strcmp(s,'Kinetic Cycles'),Araw));
numKinCyc = Eraw{ixKin};

%find 'Polarization'
iPR=find(cellfun(@(s) strcmp(s,'Polarization'),Araw));

%get all header row names
[headnum, headtxt, headraw] = xlsread(INFILE, SHEETS{snum},...
    [num2str(iPR+1) ':' num2str(iPR+1)]);
numCols = numel(headraw);
lastRow = iPR+1+numKinCyc;

%must have a header row with well names
[num,txt,raw] = xlsread(INFILE,SHEETS{snum},...
    ['A' num2str(iPR+1) ':' alphaName(numCols) num2str(lastRow)]);

%%
[numrows,numcols]=size(raw);
for col=1:numcols
    for row=1:numrows
        if strcmp(raw{row,col},'INVALID') || strcmp(raw{row,col},'Blank')
            disp(['invalid data in ' num2str(row) ', ' num2str(col)]);
            raw{row,col} = nan;
        end
    end  
end

ds=[];
ds = cell2dataset(raw);
ds.Properties.VarNames{1} = 'hours';
ds.Properties.VarNames{2} = 'time';
ds.Properties.VarNames{3} = 'temp';
ds.time = datenum((ds.time ./ 3600)) + 24*fserial{snum} - 24*fserial{1}; %time relative to first sheet
ds(isnan(ds.time),:) = [];

%need to match up column names when stacking
%they're stored in ds.Properties.VarNames
[nr,nc] = size(ds);
[nrFin,ncFin] = size(dsFinal);
for v=1:numel(ds.Properties.VarNames)
    ix=[];
    ix=find(cellfun(@(s) strcmp(s,ds.Properties.VarNames{v}),...
        dsFinal.Properties.VarNames));
    if ~isempty(ix) %found a match
        dsFinal(nrFin+(1:nr),ix) = ds(1:nr,v);
    else
        dsFinal(nrFin+(1:nr),end+1) = ds(1:nr,v);
        %if never saw this variable before, set everything above to NaN
        dsFinal(1:nrFin,end) = dataset(nan(nrFin,1)); 
        dsFinal.Properties.VarNames{end} = ds.Properties.VarNames{v};
    end
end

%deal with wells that were read in previous sheets, but not this one
for v=1:numel(dsFinal.Properties.VarNames)
    %if all entries for this column are zero in this dataset, then didn't
    %measure this well = set it to NaNs
    if sum(double(dsFinal(nrFin+(1:nr),v)) == zeros(size(nr,1))) == nr
        dsFinal(nrFin+(1:nr),v) = dataset(nan(size(nr,1)));
    end
end
%dsFinal = [dsFinal; ds];

end
end
