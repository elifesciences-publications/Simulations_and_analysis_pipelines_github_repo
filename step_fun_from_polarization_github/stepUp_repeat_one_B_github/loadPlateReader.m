function [AllRxns, AllWells, ActTime, PolData] = loadPlateReader(INFILE,INSHEETS)
%2017-02-13, EL: Load plate reader data. Read in Excel data from Tecan and
%convert it to a format used for fitting: one column for rxn names, one
%column for corresponding plate reader well names, one column for
%measurement times, another column for polarization data. Analogous to the
%way %P KaiC data was treated.
%Inputs: 'INFILE.xlsx', {'INSHEETS'}. Outputs: {AllRxns}, {AllWells},
%[ActTime], [polData].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: THE LABELS AND INDICES SPECIFIC FOR EACH EXPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

abc='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

ds = joinSheets(INFILE,INSHEETS);

%make corresponding labels for rxns and wells
for i=1:14
    sUrxns{i} = ['sU' num2str(i)];
    sUwells{i} = [abc(i+1) '' num2str(14)];
    
    sDrxns{i} = ['sD' num2str(i)];
    sDwells{i} = [abc(i+1) '' num2str(15)];
end

%the control rxns are not in the same columns
sDwells{14} = 'N13';
sUwells{14} = 'I13';

%unique names
rxns = {sUrxns{:},sDrxns{:}};
wells = {sUwells{:},sDwells{:}};

for r=1:numel(rxns)
   npts=numel(double(ds(:,wells{r})));
   AllRxns((r-1)*npts+(1:npts),1) = rxns(r); %must assign a cell array
   AllWells((r-1)*npts+(1:npts),1) = wells(r);
   ActTime((r-1)*npts+(1:npts),1) = double(ds.time(1:npts));
   PolData((r-1)*npts+(1:npts),1) = double(ds(1:npts,wells{r}));
end

end

