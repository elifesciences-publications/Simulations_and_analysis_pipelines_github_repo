%2017-02-27, EL: modify for second expt
%2017-02-13, EL: modify to work with plate reader data.

%2016-09-04, EL: generate nonparametric boostrap parameters for stepup/down
%experiment using the bestFitNonParBoot.m function. For every set of
%parameters, this function generates stepUp, stepDown fns as well as their 
%linearized versions. Saves parameters and step functions to timestamped files.

close all;
clear all;

INDIR = ['.'];
cd(INDIR);

%% enter parameter values
NUMBOOT = 1; %num. bootstraps
SAVENAME = 'repeat2_Mar30_test'; %how to name saved files?
TOSAVESTEPFUNS = 0; %save step functions?
TOSAVEPARAMS = 0; %save best fit sinusoid parameter sets?

%% fit all of the polarization data and return step functions. 
%use TORESAMPLE=0 in bestFitNonParBoot.m
for n=1:NUMBOOT
    [hiATP_param(n,:), loATP_param(n,:), STEPOUT] = bestFitNonParBoot;
    %STEPOUT = {stepUpPhase, stepUpPhaseShift, T_hiATP, sULineX, sULineY, sUlinfit, upBreak, ...
    %stepDownPhase, stepDownPhaseShift, T_loATP, sDLineX, sDLineY, sDlinfit, downBreak, resampleInd};

    [   stepUpPhase(n,:), stepUpPhaseShift(n,:), ...
        T_hiATP(n), sULineX(n,:), sULineY(n,:), sUlinfit(n,:), upBreak(n,:), ...
        stepDownPhase(n,:), stepDownPhaseShift(n,:), ...
        T_loATP(n), sDLineX(n,:), sDLineY(n,:), sDlinfit(n,:), downBreak(n,:),...
        resampleInd(n,:)] ...
        = deal(STEPOUT{:});
   
    disp(['Generating boot no. ' num2str(n)]);
end

%% compute L and D fns from these params, save into a table, save params

if TOSAVESTEPFUNS == 1
    up.phase = stepUpPhase;
    up.phaseShift = stepUpPhaseShift;
    up.breakpt = upBreak;
    up.resampleInd = resampleInd;
    down.phase = stepDownPhase;
    down.phaseShift = stepDownPhaseShift;
    down.breakpt = downBreak;
    down.resampleInd = resampleInd;
    
    %save linearized funs
    up.linPhase = sULineX;
    up.linPhaseShift = sULineY;
    up.linFit = sUlinfit;
    down.linPhase = sDLineX;
    down.linPhaseShift = sDLineY;
    down.linFit = sDlinfit;
    
    save([getDate('yyyy-mm-dd') '_stepFun_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'up','down','T_hiATP','T_loATP');
end

if TOSAVEPARAMS == 1  
    save([getDate('yyyy-mm-dd') '_params_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'hiATP_param','loATP_param','resampleInd');
end


