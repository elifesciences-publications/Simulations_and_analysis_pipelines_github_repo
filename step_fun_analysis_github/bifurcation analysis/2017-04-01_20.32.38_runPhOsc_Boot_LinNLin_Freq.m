%2017-03-15, EL: modify to do a simple bifurcation analysis as a function
%of driving frequency. Instead of driving at different dutyfractions, drive
%at different frequencies.
%2017-03-14, EL: modify for mix and match
%2017-02-27, EL: modify for merged step funs, Feb. 27
%2017-02-14, EL: modified for plate reader experiment
%2016-10-31, EL: run simulations of PO driven by linearized and full step
%functions. Save data. Do not make any plots here. 
%2016-08-22, EL: modify to make plots in phase coordinates for manuscript
%2016-08-16, EL

clear all;
close all;
clc;
INDIR = ['/home/leypunskiy/freqSweep/analysis'];
cd(INDIR);

addpath('/home/leypunskiy/MATLAB');
addpath('/home/leypunskiy/MATLAB/export_fig');

TOSAVE_DATA = 1;

%% simulation parameters
TDRIVE=[6:0.0025:48];
numcyc = [1:(100*10)]; %should be a vector, not just max value!

DUTYFRAC=0.5; %always balanced duty cycle
TEND = (max(numcyc)+0.1)*TDRIVE;

%starting phases for each dutyfrac
startphase = 0*ones(size(TDRIVE)); %[0];

%%

%how many bootstraps to evaluate?
boot = [1]; 

%load bootstrapped step funs
% NOTE! must be in 2*pi units
% INDIR=['/Users/E/Documents/Advisers/Rust/Photoperiod/'... 
%     '2017-03-14_mixMatch_LD/analysis/'];
INFILE='2017-04-01_mergedStepFuns_11.35.34.mat';
load([INDIR '/' INFILE]);

%pre-allocate save variables
peakT = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
dawnphase = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
duskphase = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
duskphaseshift = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
dawnphaseshift = zeros(numel(TDRIVE), numel(numcyc),numel(boot));

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
    
%     %USE linearized versions
%     STEPUP.phase = up.linPhase(1,:);
%     STEPUP.phaseShift = up.linPhaseShift(1,:);
%     STEPDOWN.phase = down.linPhase(1,:);
%     STEPDOWN.phaseShift = down.linPhaseShift(1,:);

    STEPUP.linPhase = up.linPhase(1,:);
    STEPUP.linPhaseShift = up.linPhaseShift(1,:);
    STEPDOWN.linPhase = down.linPhase(1,:);
    STEPDOWN.linPhaseShift = down.linPhaseShift(1,:);
   
    %time this
    tictime=tic;
    clc; %clear screen after every simulation
    disp(['Starting boot ' num2str(b) '...']);
    for fr=1:numel(TDRIVE)
            tFr = tic;
            %catch errors in integration due to kinks
            try
                disp(['fr = ' num2str(TDRIVE(fr))]);
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot_Fast_Two(TDRIVE(fr), TLIGHT, TDARK,...
                    STEPUP, STEPDOWN, DUTYFRAC, TEND(fr), ...
                    numcyc, startphase(fr));
                
                peakT(fr, 1:numel(numcyc),b) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
                dawnphase(fr,1:numel(numcyc),b) = DAWNPHASE;
                duskphase(fr,1:numel(numcyc),b) = DUSKPHASE;
                duskphaseshift(fr,1:numel(numcyc),b) = DUSKPHASESHIFT;
                dawnphaseshift(fr,1:numel(numcyc),b) = DAWNPHASESHIFT;
            catch err
                disp(['boot ' num2str(b) ...
                    ', fr ' num2str(TDRIVE(fr)) ...
                     ' simulation failed']);
                msgText = getReport(err);
                warning(['ErrorMsg:' msgText]);
                peakT(fr, 1:numel(numcyc),b) = nan;
                dawnphase(fr,1:numel(numcyc),b) = nan;
                duskphase(fr,1:numel(numcyc),b) = nan;
                duskphaseshift(fr,1:numel(numcyc),b) = nan;
                dawnphaseshift(fr,1:numel(numcyc),b) = nan;
            end
            
            disp(['finished in ' num2str(toc(tFr)/60) ' min.' char(10)]);
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% save simulation results
if TOSAVE_DATA == 1
    OUTDIR='.';
    SAVENAME=['freqSweep' '_stepFunMixMatch_' 'NonLin_' num2str(max(boot))];
    NPBoot.peakT = peakT;
    NPBoot.dawnphase = dawnphase;
    NPBoot.dawnphaseshift = dawnphaseshift;
    NPBoot.duskphase = duskphase;
    NPBoot.duskphaseshift = duskphaseshift;
    NPBoot.numcyc = numcyc;
    NPBoot.dutyfrac = DUTYFRAC;
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

break;

%pre-allocate save variables for step funs
peakT = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
dawnphase = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
duskphase = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
duskphaseshift = zeros(numel(TDRIVE), numel(numcyc),numel(boot));
dawnphaseshift = zeros(numel(TDRIVE), numel(numcyc),numel(boot));


%% repeat for linearized funs: run simulations with different L and D funs
for b=boot
    up = up_mix{b};
    down = down_mix{b};
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    
    %using slightly different syntax than in previous files -- here, choose
    %up and down funs from the cell arrays above, but each contains only
    %one line of phases/phaseShifts, so replace all row indices with 1.
%     STEPUP.phase = up.phase(1,:);
%     STEPUP.phaseShift = up.phaseShift(1,:);
%     STEPDOWN.phase = down.phase(1,:);
%     STEPDOWN.phaseShift = down.phaseShift(1,:);
    
    %USE linearized versions
    STEPUP.phase = up.linPhase(1,:);
    STEPUP.phaseShift = up.linPhaseShift(1,:);
    STEPDOWN.phase = down.linPhase(1,:);
    STEPDOWN.phaseShift = down.linPhaseShift(1,:);
    
    STEPUP.linPhase = up.linPhase(1,:);
    STEPUP.linPhaseShift = up.linPhaseShift(1,:);
    STEPDOWN.linPhase = down.linPhase(1,:);
    STEPDOWN.linPhaseShift = down.linPhaseShift(1,:);
   
    %time this
    tictime=tic;
    disp(['Starting linearized boot ' num2str(b) '...']);
    for fr=1:numel(TDRIVE)
            tStart = tic;         
            %catch errors in integration due to kinks
            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot_Fast(TDRIVE(fr), TLIGHT, TDARK,...
                    STEPUP, STEPDOWN, DUTYFRAC, TEND(fr), ...
                    numcyc, startphase(fr));
                
                peakT(fr, 1:numel(numcyc),b) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
                dawnphase(fr,1:numel(numcyc),b) = DAWNPHASE;
                duskphase(fr,1:numel(numcyc),b) = DUSKPHASE;
                duskphaseshift(fr,1:numel(numcyc),b) = DUSKPHASESHIFT;
                dawnphaseshift(fr,1:numel(numcyc),b) = DAWNPHASESHIFT;
            catch err
                disp(['boot ' num2str(b) ...
                    ', fr ' num2str(TDRIVE(fr)) ...
                     ' simulation failed']);
                peakT(fr, 1:numel(numcyc),b) = nan;
                dawnphase(fr,1:numel(numcyc),b) = nan;
                duskphase(fr,1:numel(numcyc),b) = nan;
                duskphaseshift(fr,1:numel(numcyc),b) = nan;
                dawnphaseshift(fr,1:numel(numcyc),b) = nan;
            end
            tEnd = toc(tStart)/60;
            disp(['fr=' num2str(TDRIVE(fr)) ' finished in ' num2str(tEnd) ' min.']); 	
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Linearized boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% save simulation results
if TOSAVE_DATA == 1
    OUTDIR='.';
    SAVENAME=['freqSqeep' '_stepFunMixMatch_' 'Lin_' num2str(max(boot))];
    
    NPBoot.peakT = peakT;
    NPBoot.dawnphase = dawnphase;
    NPBoot.dawnphaseshift = dawnphaseshift;
    NPBoot.duskphase = duskphase;
    NPBoot.duskphaseshift = duskphaseshift;
    NPBoot.numcyc = numcyc;
    NPBoot.dutyfrac = DUTYFRAC;
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
