%2107-03-30, EL: perform a global fit to all of the data from an experiment
%measuring D(theta) step-response function. Save data. 
%
% set TORESAMPLE=0 in bestFitNonParBoot.m.
%

close all;
clear all;

%% load raw data, step times; set normalization/fit parameters
TOSAVEPARAMS = 1;
SAVENAME = 'Feb28_plateReader_FullFile_';
TOSAVESTEPFUNS = 1;

INDIR = ['.'];
cd(INDIR);

%% draw nonparametric boostraps by resampling %P KaiC data (pNet_trim vector)
NUMBOOT = 1; %num. bootstraps
for n=1:NUMBOOT
    [hiATP_param(n,:), loATP_param(n,:), STEPOUT] = bestFitNonParBoot;
    %STEPOUT = {stepUpPhase, stepUpPhaseShift, T_hiATP, sULineX, sULineY, sUlinfit, upBreak, ...
    %stepDownPhase, stepDownPhaseShift, T_loATP, sDLineX, sDLineY, sDlinfit, downBreak, resampleInd};

    [~, ~, ...
        T_hiATP(n), ~,~, ~, ~, ...
        stepDownPhase(n,:), stepDownPhaseShift(n,:), ...
        T_loATP(n), sDLineX(n,:), sDLineY(n,:), sDlinfit(n,:), downBreak(n,:),...
        resampleInd(n,:)] ...
        = deal(STEPOUT{:});
   
    disp(['Generating boot no. ' num2str(n)]);
end

%% compute L and D fns from these params, save into a table, save params
if TOSAVESTEPFUNS == 1
    down.phase = stepDownPhase;
    down.phaseShift = stepDownPhaseShift;
    down.breakpt = downBreak;
    down.resampleInd = resampleInd;
    down.linPhase = sDLineX;
    down.linPhaseShift = sDLineY;
    down.linFit = sDlinfit;
    
    save([getDate('yyyy-mm-dd') '_stepFun_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'down','T_hiATP','T_loATP');
end

if TOSAVEPARAMS == 1  
    save([getDate('yyyy-mm-dd') '_params_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'hiATP_param','loATP_param','resampleInd');
end


