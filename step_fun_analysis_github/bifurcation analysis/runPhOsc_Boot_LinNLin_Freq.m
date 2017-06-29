%2017-03-28, EL: run simulations of a phase oscillator subjected to light-dark
%cycles of different period. 

clear all;
close all;
clc;
cd(['.']);

TOSAVE_DATA = 0;

%% simulation parameters
TDRIVE=[6:1:48]; %driving periods
                 %figure in text: TDRIVE=[6:0.0025:48];
                 
%after how many cycles do you want to store phases at dawn/dusk and phase shifts?
numcyc = [1:100]; %muct be a vector, not just max value!
                  %figure in text: numcyc = [1:(100*10)]; 

DUTYFRAC=0.5; %always balanced duty cycle
TEND = (max(numcyc)+0.1)*TDRIVE;

%starting phases for each dutyfrac
startphase = 0*ones(size(TDRIVE)); %[0];

%%

%how many bootstraps to evaluate?
boot = [1]; 

%load bootstrapped step funs
% NOTE! must be in 2*pi units
INDIR=['../helper functions and shared files/'];
INFILE='2017-06-05_widefit_mergedStepFuns_11.30.54.mat';
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
            disp(char(10));
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