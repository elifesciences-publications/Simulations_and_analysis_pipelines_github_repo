
clear all;
close all;
clc;
INDIR = ['.'];
cd(INDIR);

TOSAVE_DATA = 1;

%% simulation parameters
dutyfrac=[4:1:18]/24;
numcyc = [1:20]; %should be a vector, not just max value!

TDRIVE=24.0;
TEND = (max(numcyc)+0.1)*TDRIVE;

%start from ph=0
startphase = zeros(size(dutyfrac)); %[0]; 

%%

%how many bootstraps to evaluate?
boot = [1:4]; 

%load bootstrapped step funs
% NOTE! must be in 2*pi units
INDIR=['helper functions and shared files'];
% INFILE='2017-03-17_stepFun_inKaiCPhase_19.20.23.mat';
% INFILE='2017-03-30_mergedStepFuns_11.43.53.mat';
% INFILE='2017-03-31_mergedStepFuns_21.04.38.mat';
INFILE='2017-04-01_mergedStepFuns_11.35.34.mat';
load(['../' INDIR '/' INFILE]);

%% run simulations with different L and D funs
for b=boot
    up = up_mix{b};
    down = down_mix{b};
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    
    %using slightly different syntax than in previous files -- here, choose
    %up and down funs from the cell arrays above, but each contains only
    %one line of phases/phaseShifts, so replace all row indices with 1.
    STEPUP.phase = up.phase(1,:);
    STEPUP.phaseShift = up.phaseShift(1,:);
    STEPDOWN.phase = down.phase(1,:);
    STEPDOWN.phaseShift = down.phaseShift(1,:);
   
    STEPUP.linPhase = up.linPhase(1,:);
    STEPUP.linPhaseShift = up.linPhaseShift(1,:);
    STEPDOWN.linPhase = down.linPhase(1,:);
    STEPDOWN.linPhaseShift = down.linPhaseShift(1,:);
   
    %time this
    tictime=tic;
    disp(['Starting boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
            disp(['df=' num2str(df)]);

            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot_Fast_Two(TDRIVE, TLIGHT, TDARK,...
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
    SAVENAME=['NonLin_FAST_' num2str(max(boot))];
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
    up = up_mix{b};
    down = down_mix{b};
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    
    %using slightly different syntax than in previous files -- here, choose
    %up and down funs from the cell arrays above, but each contains only
    %one line of phases/phaseShifts, so replace all row indices with 1
    
    %USE linearized versions
    STEPUP.phase = up.linPhase(1,:);
    STEPUP.phaseShift = up.linPhaseShift(1,:);
    STEPDOWN.phase = down.linPhase(1,:);
    STEPDOWN.phaseShift = down.linPhaseShift(1,:);
   
    %time this
    tictime=tic;
    disp(['Starting linearized boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot_Fast_Two(TDRIVE, TLIGHT, TDARK,...
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
    disp(['Linearized boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% save simulation results
if TOSAVE_DATA == 1
    OUTDIR='.';
    SAVENAME=['Lin_FAST_' num2str(max(boot))];
    
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
    
    saveMyData(NPBoot,SAVENAME,OUTDIR);
    
    %save this script with timestamp
    thisFileName=[mfilename('fullpath')];
    splitFileName = strsplit(thisFileName,'/');
    datedName = [getDate() '_' splitFileName{end}];
    saveFileName=[OUTDIR '/' datedName '.m'];
    currentFile=strcat(thisFileName, '.m');
    copyfile(currentFile,saveFileName);
    
end
