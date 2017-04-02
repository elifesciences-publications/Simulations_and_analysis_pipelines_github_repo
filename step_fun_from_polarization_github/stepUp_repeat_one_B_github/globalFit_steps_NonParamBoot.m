%2017-03-30, EL: generate the L step function from the 05-Feb-2017
%polarization measurements.

close all;
clear all;
INDIR = ['.'];
cd(INDIR);

%% enter input parameters here
SAVENAME = 'Feb_plateReader_Apr1analyze_'; %name to append to savefiles?
TOSAVEPARAMS = 1; %save fit parameters? (1=yes)
TOSAVESTEPFUNS = 1; %save step functions? (1=yes)

%% fit all polarization data and get stepUp/stepDown funs
NUMBOOT = 1; %num. bootstraps to draw (use 1 here and set TORESAMPLE=0 in bestFitNonParBoot.m)

for n=1:NUMBOOT
    [hiATP_param(n,:), ~, STEPOUT] = bestFitNonParBoot;
    %STEPOUT = {stepUpPhase, stepUpPhaseShift, T_hiATP, sULineX, sULineY, sUlinfit, upBreak, ...
    %stepDownPhase, stepDownPhaseShift, T_loATP, sDLineX, sDLineY, sDlinfit, downBreak, resampleInd};

    [stepUpPhase(n,:), stepUpPhaseShift(n,:), ...
        T_hiATP(n), sULineX(n,:), sULineY(n,:), sUlinfit(n,:), upBreak(n,:), ...
        ~, ~, ...
        T_loATP(n), ~, ~, ~, ~,...
        resampleInd(n,:)] ...
        = deal(STEPOUT{:});
   
    disp(['Generating boot no. ' num2str(n)]);
end

%% compute L fn from these params, save into a table, save params

if TOSAVESTEPFUNS == 1
    up.phase = stepUpPhase;
    up.phaseShift = stepUpPhaseShift;
    up.breakpt = upBreak;
    up.resampleInd = resampleInd;
       
    %save linearized funs
    up.linPhase = sULineX;
    up.linPhaseShift = sULineY;
    up.linFit = sUlinfit;
        
    save([getDate('yyyy-mm-dd') '_stepFun_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'up', 'T_hiATP','T_loATP');
end

if TOSAVEPARAMS == 1  
    save([getDate('yyyy-mm-dd') '_params_' SAVENAME num2str(NUMBOOT) '_' getDate() '.mat'], ...
        'hiATP_param','resampleInd');
end


