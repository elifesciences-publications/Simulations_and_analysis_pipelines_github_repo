%2017-03-29, EL: run phase oscillator simulations using bootstrapped
%step-response functions based on KaiC phosphorylation data.
%
% Load bootstrapped L, D, Llin, Dlin, TDark, TLight values from INFILE. 
% You can use the globalFit_steps_NonParamBoot.m script to generate such a
% file. The file with 1000 bootstrapped parameter sets used in Fig.
% 4-figSup2 simulations is located in the saved_data/ folder: 
% INFILE='2016-10-31_stepFun_Lin-NLin-Oct_1000_2016-10-31_20.45.44.mat'.
%
% Set simulation parameters below. 

clear all;
close all;
clc;
cd(['.']);

TOSAVE_DATA = 1; %save timestamped simulation data and timestamped script?

%% simulation parameters
TDRIVE=24.0; %driving period (hrs)
numcyc = 1:10; %after which light-dark cycles to store dawnphase, duskphase, etc.
TEND = max(numcyc)*TDRIVE+10; 
dutyfrac=[4:0.25:18]/24; %what day lengths to simulate? (as frac of TDRIVE)

%compute startphase based on entrained conditions observed in vitro
inVitroDawnPhase =  [0.2525    0.1832    0.1486    0.1038    0.0676   -0.0080];
inVitroDuty = [0.2500    0.3333    0.4167    0.5000    0.5833    0.7500];
startphase = interp1(inVitroDuty, inVitroDawnPhase, dutyfrac,...
'linear','extrap');

%which bootstrapped L and D functions to use in simulations?
boot = [1:10]; 

%load bootstrapped step funs
% NOTE! must be in 2*pi units
INDIR=['saved_data'];
INFILE='2017-03-29_stepFun_Lin-NLin-_20_22.22.18.mat';
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
                     'simulation failed']);
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
    SAVENAME=['NonParamBoot_' num2str(max(boot))];
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
    TLIGHT=T_hiATP(b);
    TDARK=T_loATP(b);
    
    %USE linearized versions
    STEPUP.phase = up.linPhase(b,:);
    STEPUP.phaseShift = up.linPhaseShift(b,:);
    STEPDOWN.phase = down.linPhase(b,:);
    STEPDOWN.phaseShift = down.linPhaseShift(b,:);
  
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
                     'simulation failed']);
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
    SAVENAME=['NonParamBoot_Lin_' num2str(max(boot))];
    
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
