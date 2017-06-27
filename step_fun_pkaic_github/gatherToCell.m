function [rxnx, rxny, rxnx_tofit, rxny_tofit, rxnstep, rxnlabel, ...
    rxnmu, rxnsigma] = ...
    gatherToCell(rxnset, steptime, afterstep, AllRxns, ActTime, pNet, TONORM)
%collect stepUp/stepDown datasets into cell arrays for fitting
for r=1:numel(rxnset)
    ind = cellfun(@(x) strcmp(x,rxnset{r}), AllRxns);
    t=ActTime(ind);
    y=pNet(ind);
     
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
    
    % keep all of rxn if it's a control
    if strcmp(rxnset{r},'sD10') || strcmp(rxnset{r},'sU10')
        rxnx_tofit{r} = t;
        rxny_tofit{r} = y;
        rxnstep{r} = nan;
        rxnlabel{r} = 'control';
    else 
        % begin fit afterstep hrs after step 
        goodt = (t >= steptime(r)+afterstep); %fixed to >= from >; should include step pt if afterstep=0
        rxnx_tofit{r} = t(goodt);
        rxny_tofit{r} = y(goodt);
        rxnstep{r} = steptime(r);
        
        if strcmp(rxnset{r}(1:2), 'sU')
            rxnlabel{r} = ['stepUp at t=' num2str(steptime(r))];
        elseif strcmp(rxnset{r}(1:2),'sD')
            rxnlabel{r} = ['stepDown at t=' num2str(steptime(r))];
        end
    end
    
end

% disp(['in gatherToCell:']);
% disp(['rxnx_tofit: ' rxnx_tofit]);
% disp(['rxny_tofit: ' rxny_tofit]);
% disp(['rxnx: ' rxnx]);
% disp(['rxny: ' rxny]);


end