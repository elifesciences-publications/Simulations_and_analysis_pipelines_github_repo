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
for i=1:16
    if i < 9
        sDrxns{i} = ['sD' num2str(i)];
        sDwells{i} = ['F' num2str(i+13)];
        sUrxns{i} = ['sU' num2str(i)];
        sUwells{i} = ['H' num2str(i+13)];
    elseif i >= 9 && i <=16
        sDrxns{i} = ['sD' num2str(i)];
        sDwells{i} = ['G' num2str(i+5)];
        sUrxns{i} = ['sU' num2str(i)];
        sUwells{i} = ['I' num2str(i+5)];
        if i == 16
            sUwells{i} = ['I22']; 
        end
    end
end


%unique names
rxns = {sDrxns{:},sUrxns{:}};
wells = {sDwells{:},sUwells{:}};

for r=1:numel(rxns)
   npts=numel(double(ds(:,wells{r})));
   AllRxns((r-1)*npts+(1:npts),1) = rxns(r); %must assign a cell array
   AllWells((r-1)*npts+(1:npts),1) = wells(r);
   ActTime((r-1)*npts+(1:npts),1) = double(ds.time(1:npts));
   PolData((r-1)*npts+(1:npts),1) = double(ds(1:npts,wells{r}));
end

end

