%2017-05-29, EL: revise to account for period differences in hi vs lo ATP
%2016-09-04, EL: generate nonparametric boostrap parameters for stepup/down
%experiment using the bestFitNonParBoot.m function. For every set of
%parameters, this function generates stepUp, stepDown fns as well as their 
%linearized versions. Saves parameters and step functions to timestamped files.

close all;
clear all;

%% input parameters
INDIR = ['.'];
cd(INDIR);
TOSAVEPARAMS = 1; %Save fit parameters for every bootstrap? (1=yes)
TOSAVESTEPFUNS = 1; %Save bootstrapped step functions? (1=yes)
SAVENAME = 'JuneFit_updFitPipe_wideLinFit';

%% draw nonparametric boostraps by resampling %P KaiC data 
NUMBOOT = 1000; %how many bootstraps to generate?

for n=1:NUMBOOT
    [hiATP_param(n,:), loATP_param(n,:), STEPOUT] = bestFitNonParBoot;
    %[hiATP_param(n,:), loATP_param(n,:), STEPOUT] = bestFitNonParBoot;
    %STEPOUT = {stepUpPhase, stepUpPhaseShift, T_hiATP, sULineX, sULineY, sUlinfit, upBreak, ...
    %stepDownPhase, stepDownPhaseShift, T_loATP, sDLineX, sDLineY, sDlinfit, downBreak, resampleInd};

    [stepUpPhase(n,:), stepUpPhaseShift(n,:), ...
        T_hiATP(n), sULineX(n,:), sULineY(n,:), sUlinfit(n,:), upBreak(n,:) ...
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
    
    save([getDate('yyyy-mm-dd') '_stepFun_' SAVENAME num2str(NUMBOOT) ...
        '_' getDate('HH.MM.SS') '.mat'], ...
        'up', 'down','T_hiATP','T_loATP');
end


if TOSAVEPARAMS == 1  
    save([getDate('yyyy-mm-dd') '_params_' SAVENAME num2str(NUMBOOT) ...
        '_' getDate('HH.MM.SS') '.mat'], ...
        'hiATP_param','loATP_param','resampleInd');
end


