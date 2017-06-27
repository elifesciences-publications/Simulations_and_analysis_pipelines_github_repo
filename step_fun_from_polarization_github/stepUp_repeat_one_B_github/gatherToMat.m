function [ActTime_trim, PolData_trim, AllRxns_trim, AllWells_trim] = ...
    gatherToMat(rxnset, delayTime, AllRxns, AllWells, ActTime, PolData, TONORM)
%2017-02-13, EL: modified for plate reader data that's in LL effectively.
%So just toss data before delayTime hours after starting the measurement.

%collect stepUp/stepDown datasets into an array that contains only
%datapoints used for fitting. output of this function will be used to draw
%resampled datasets. 
%probably shouldn't normalize raw data, but only data used for fits!


ActTime_trim = [];
PolData_trim = [];
AllRxns_trim = [];
AllWells_trim = [];

for r=1:numel(rxnset)
    ind = cellfun(@(x) strcmp(x,rxnset{r}), AllRxns);
    t=ActTime(ind);
    y=PolData(ind);
    names=AllRxns(ind);
    wells=AllWells(ind);
    
%     disp(['r=' num2str(r)]);
%     disp('actTime trim: ');
%     disp(ActTime_trim);
%     disp('pNet trim: ');
%     disp(pNet_trim);
%     disp('AllRxns_trim: ');
%     disp(AllRxns_trim);
%     disp('Names:');
%     disp(names);
    
    % mean0, unit variance normalization
    %probably shouldn't normalize raw data, but only data used for fits!
    %MOVE THIS TO END?!
    if TONORM==1
        y = 100*(y - mean(y))./std(y);
    end
    
    % keep all of rxn if it's a control
    % begin fit afterstep hrs after step
    goodt = (t > delayTime);
    savet = t(goodt);
    savey = y(goodt);
    savenames = names(goodt);
    savewells = wells(goodt);
    
    if isempty(ActTime_trim) & isempty(PolData_trim) & isempty(AllRxns_trim)
        ActTime_trim = savet;
        PolData_trim = savey;
        AllRxns_trim = savenames;
        AllWells_trim = savewells;
    else
        ActTime_trim = [ActTime_trim; savet];
        PolData_trim = [PolData_trim; savey];
        AllRxns_trim = {AllRxns_trim{:}, savenames{:}};
        AllWells_trim = {AllWells_trim{:}, savewells{:}};
    end  
    
end
end