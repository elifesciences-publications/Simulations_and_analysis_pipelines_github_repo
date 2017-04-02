function [rxnx, rxny, rxnx_tofit, rxny_tofit, rxnstep, rxnlabel, ...
    rxnmu, rxnsigma] = ...
    gatherToCell(rxnset, steptime, delayTime, AllRxns, ActTime, PolData, TONORM)
%2017-02-13, EL: modified for plate reader data
%collect stepUp/stepDown datasets into cell arrays for fitting
for r=1:numel(rxnset)
    y=[];
    t=[];
    ind=[];
    
    ind = cellfun(@(x) strcmp(x,rxnset{r}), AllRxns);
    t=ActTime(ind);
    y=PolData(ind);
    
    bady = isnan(y);
    y=y(~bady);
    t=t(~bady);
     
    % mean0, unit variance normalization
    if TONORM==1
        rxnmu{r} = mean(y);
        rxnsigma{r} = std(y);
        y = 100*(y - mean(y))./std(y);
    else
        rxnmu{r} = nan;
        rxnsigma{r} = nan;
        %y = 100*y; %still put in %P KaiC
    end
    
    rxnx{r} = t; %raw
    rxny{r} = y; %raw
    
    % begin fit delayTime hrs after putting on plate reader
    goodt = (t >= delayTime); %fixed to >= from >; should include step pt if afterstep=0
    rxnx_tofit{r} = t(goodt);
    rxny_tofit{r} = y(goodt);
    rxnstep{r} = steptime(r);
    
    % make labels
    if strcmp(rxnset{r},'sD8') || strcmp(rxnset{r},'sD16')
        rxnstep{r} = nan;
        rxnlabel{r} = 'control';
    else  
        if strcmp(rxnset{r}(1:2), 'sU')
            rxnlabel{r} = ['stepUp at t=' num2str(steptime(r),'%2.1f')];
        elseif strcmp(rxnset{r}(1:2),'sD')
            rxnlabel{r} = ['stepDown at t=' num2str(steptime(r),'%2.1f')];
        end
    end
    
end

% disp(['in gatherToCell:']);
% disp(['rxnx_tofit: ' rxnx_tofit]);
% disp(['rxny_tofit: ' rxny_tofit]);
% disp(['rxnx: ' rxnx]);
% disp(['rxny: ' rxny]);


end