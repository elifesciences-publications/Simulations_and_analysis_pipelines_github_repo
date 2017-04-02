
clear all;
close all;
clc;
cd(['.']);

TOSAVE_DATA = 1;

%% simulation parameters
dutyfrac=[6:1:18]/24;
numcyc = [1:50]; %should be a vector, not just max value!

TDRIVE=24.0;
TEND = (max(numcyc)+0.1)*TDRIVE;

startphase = zeros(size(dutyfrac)); 

%%

%how many bootstraps to evaluate?
boot = [1]; 

%load bootstrapped step funs
% NOTE! must be in 2*pi units
INDIR=['saved_data/'];
INFILE='2017-03-30_stepFun_repeat2_Mar30_test1_2017-03-30_09.53.38.mat';
load([INDIR '/' INFILE]);

%% run simulations with different L and D funs
for b=boot
    TLIGHT=T_hiATP(b);
    TDARK=T_loATP(b);
    STEPUP.phase = up.phase(b,:);
    STEPUP.phaseShift = up.phaseShift(b,:);
    STEPDOWN.phase = down.phase(b,:);
    STEPDOWN.phaseShift = down.phaseShift(b,:);
    
    STEPUP.linPhase = up.linPhase(b,:);
    STEPUP.linPhaseShift = up.linPhaseShift(b,:);
    STEPDOWN.linPhase = down.linPhase(b,:);
    STEPDOWN.linPhaseShift = down.linPhaseShift(b,:);
   
    %time this
    tictime=tic;
    disp(['Starting boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
        
            %catch errors in integration due to kinks
            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot(TDRIVE, TLIGHT, TDARK,...
                    STEPUP, STEPDOWN, dutyfrac(df), TEND, ...
                    numcyc, startphase(df));
                
                peakT(df, 1:numel(numcyc),b) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
                dawnphase(df,1:numel(numcyc),b) = DAWNPHASE;
                duskphase(df,1:numel(numcyc),b) = DUSKPHASE;
                duskphaseshift(df,1:numel(numcyc),b) = DUSKPHASESHIFT;
                dawnphaseshift(df,1:numel(numcyc),b) = DAWNPHASESHIFT;
            catch err
                disp(['boot ' num2str(b) ...
                    ', df ' num2str(dutyfrac(df)) ...
                     ' simulation failed']);
                peakT(df, 1:numel(numcyc),b) = nan;
                dawnphase(df,1:numel(numcyc),b) = nan;
                duskphase(df,1:numel(numcyc),b) = nan;
                duskphaseshift(df,1:numel(numcyc),b) = nan;
                dawnphaseshift(df,1:numel(numcyc),b) = nan;
            end
        
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% save simulation results
if TOSAVE_DATA == 1
    OUTDIR='.';
    SAVENAME=['NonParamBoot_I22_HiBrPt_' num2str(max(boot))];
    NPBoot.peakT = peakT;
    NPBoot.dawnphase = dawnphase;
    NPBoot.dawnphaseshift = dawnphaseshift;
    NPBoot.duskphase = duskphase;
    NPBoot.duskphaseshift = duskphaseshift;
    NPBoot.numcyc = numcyc;
    NPBoot.dutyfrac = dutyfrac;
    NPBoot.boot = boot;
    NPBoot.startphase = startphase;
    NPBoot.TDRIVE = TDRIVE;
    NPBoot.TEND = TEND;
    NPBoot.numcyc = numcyc;
    
    saveMyData(NPBoot,SAVENAME,OUTDIR);
        
    %save this script with timestamp
    thisFileName=[mfilename('fullpath')];
    splitFileName = strsplit(thisFileName,'/');
    datedName = [getDate() '_' splitFileName{end}];
    saveFileName=[OUTDIR '/' datedName '.m'];
    currentFile=strcat(thisFileName, '.m');
    copyfile(currentFile,saveFileName);
    
    clear peakT dawnphase duskphase duskphaseshift dawnphaseshift NPBoot;
    
end

%% repeat for linearized funs: run simulations with different L and D funs
for b=boot
    TLIGHT=T_hiATP(b);
    TDARK=T_loATP(b);
    
    %USE linearized versions
    STEPUP.phase = up.linPhase(b,:);
    STEPUP.phaseShift = up.linPhaseShift(b,:);
    STEPDOWN.phase = down.linPhase(b,:);
    STEPDOWN.phaseShift = down.linPhaseShift(b,:);

    STEPUP.linPhase = up.linPhase(b,:);
    STEPUP.linPhaseShift = up.linPhaseShift(b,:);
    STEPDOWN.linPhase = down.linPhase(b,:);
    STEPDOWN.linPhaseShift = down.linPhaseShift(b,:);
   
    %time this
    tictime=tic;
    disp(['Starting linearized boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
         
            %catch errors in integration due to kinks
            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot(TDRIVE, TLIGHT, TDARK,...
                    STEPUP, STEPDOWN, dutyfrac(df), TEND, ...
                    numcyc, startphase(df));
                
                peakT(df, 1:numel(numcyc),b) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
                dawnphase(df,1:numel(numcyc),b) = DAWNPHASE;
                duskphase(df,1:numel(numcyc),b) = DUSKPHASE;
                duskphaseshift(df,1:numel(numcyc),b) = DUSKPHASESHIFT;
                dawnphaseshift(df,1:numel(numcyc),b) = DAWNPHASESHIFT;
            catch err
                disp(['boot ' num2str(b) ...
                    ', df ' num2str(dutyfrac(df)) ...
                     ' simulation failed']);
                peakT(df, 1:numel(numcyc),b) = nan;
                dawnphase(df,1:numel(numcyc),b) = nan;
                duskphase(df,1:numel(numcyc),b) = nan;
                duskphaseshift(df,1:numel(numcyc),b) = nan;
                dawnphaseshift(df,1:numel(numcyc),b) = nan;
            end
        
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% save simulation results
if TOSAVE_DATA == 1
    OUTDIR='.';
    SAVENAME=['NonParamBoot_Lin_I22_HiBrPt_' num2str(max(boot))];
    
    NPBoot.peakT = peakT;
    NPBoot.dawnphase = dawnphase;
    NPBoot.dawnphaseshift = dawnphaseshift;
    NPBoot.duskphase = duskphase;
    NPBoot.duskphaseshift = duskphaseshift;
    NPBoot.numcyc = numcyc;
    NPBoot.dutyfrac = dutyfrac;
    NPBoot.boot = boot;
    NPBoot.startphase = startphase;
    NPBoot.TDRIVE = TDRIVE;
    NPBoot.TEND = TEND;
    NPBoot.numcyc = numcyc;
    
    saveMyData(NPBoot,SAVENAME,OUTDIR);
    
    %save this script with timestamp
    thisFileName=[mfilename('fullpath')];
    splitFileName = strsplit(thisFileName,'/');
    datedName = [getDate() '_' splitFileName{end}];
    saveFileName=[OUTDIR '/' datedName '.m'];
    currentFile=strcat(thisFileName, '.m');
    copyfile(currentFile,saveFileName);
    
end
