function [ActTime_trim, pNet_trim, AllRxns_trim] = ...
    gatherToMat(rxnset, steptime, afterstep, AllRxns, ActTime, pNet, TONORM)
%collect stepUp/stepDown datasets into an array that contains only
%datapoints used for fitting. output of this function will be used to draw
%resampled datasets. 
%probably shouldn't normalize raw data, but only data used for fits!

ActTime_trim = [];
pNet_trim = [];
AllRxns_trim = [];

for r=1:numel(rxnset)
    ind = cellfun(@(x) strcmp(x,rxnset{r}), AllRxns);
    t=ActTime(ind);
    y=pNet(ind);
    names=AllRxns(ind);
    
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
    if strcmp(rxnset{r},'sD10') || strcmp(rxnset{r},'sU10')
        if isempty(ActTime_trim) & isempty(pNet_trim) & isempty(AllRxns_trim)
           ActTime_trim = t;
           pNet_trim = y;
           AllRxns_trim = names;
        else
            ActTime_trim = [ActTime_trim; t];
            pNet_trim = [pNet_trim; y];
            AllRxns_trim = {AllRxns_trim{:}, names{:}};
        end
    else 
        % begin fit afterstep hrs after step 
        goodt = (t > steptime(r)+afterstep);
        savet = t(goodt);
        savey = y(goodt);
        savenames = names(goodt);
        
        if isempty(ActTime_trim) & isempty(pNet_trim) & isempty(AllRxns_trim)
           ActTime_trim = savet;
           pNet_trim = savey;
           AllRxns_trim = savenames;
        else
            ActTime_trim = [ActTime_trim; savet];
            pNet_trim = [pNet_trim; savey];
            AllRxns_trim = {AllRxns_trim{:}, savenames{:}};
        end
        
    end
    
end
end