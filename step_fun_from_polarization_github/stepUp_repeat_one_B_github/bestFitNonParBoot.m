%2017-05-29, EL: revised calculations
%2017-02-13, EL: modify for plate reader data
%2016-09-04, EL: perform one non-parametric bootstrap resampling of
%stepup/down data and fit. 
function [PARSET_hiATP, PARSET_loATP, STEPOUT] = ...
    bestFitNonParBoot
%Return a set of best fit parameters for a set of resampled loATP and hiATP
%datasets, as well as step functions (and their linearizations) that result
%from this best fit parameter set.

%% load raw data, step times; set normalization/fit parameters
INDIR = ['saved_data'];
INFILE = [INDIR '/2017-02-05_stepUpDown_plateReader_forMatlab.xlsx'];
INSHEETS = {'Sheet6','Sheet7'};

TORESAMPLE = 0; %if 0, don't bootstrap -- use entire dataset
TOTEST = 1; %test goodness of fits by plotting? (1=yes)
TOPLOT_STEPFUNS = 1; %plot step functions? (1=yes)
TOEXPORT_fBootStep = 1; %export figures? (1=yes)

% script to load plate reader data
[AllRxns, AllWells, ActTime, PolData] = loadPlateReader(INFILE,INSHEETS);

%'stepTimes' sheet has stepTimes relative to t=0 on plate reader
[stepTimes,~,~] = xlsread(INFILE,'stepTimes','C2:C14');
stepUpTime = [stepTimes; nan]; 
stepDownTime = [stepTimes; nan];

afterstep = 2; %how much time to wait after step to start fitting
TONORM = 1; %normalize to mean0,std100 or to keep as raw data for fits/plots
TONORMTRIM = 0; %don't normalize while trimming
%%
allRxnTypes = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9',...
    'sU10','sU11','sU12','sU13','sU14',...
               'sD14'};
         
[ActTime_trim, PolData_trim, AllRxns_trim, AllWells_trim] = ...
    gatherToMat(allRxnTypes, 2, AllRxns, AllWells,ActTime, PolData, TONORMTRIM);
AllRxns_trim = AllRxns_trim'; %must transpose
AllWells_trim = AllWells_trim';

%% draw a random set of datapoints with replacement
numpts = numel(PolData_trim); %total num. pts

%%%%%%%%%%%
numnew = 1; %MUST BE 1 here. generate one fit set at a time
%%%%%%%%%%%

%random draw
npBootInd = randi(numpts, numpts, numnew);
npBootInd = sort(npBootInd);

%allow no resampling
if TORESAMPLE == 0
    npBootInd = 1:numpts;
end

%draw a single set of data
AllRxns_trim = AllRxns_trim(npBootInd);
PolData_trim = PolData_trim(npBootInd);
ActTime_trim = ActTime_trim(npBootInd);

%% gather new datapoints into cell arrays
%put all high-ATP rxns into a single cell array
hiATPset = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9','sU10',...
    'sU11','sU12','sU13','sD14'};
hiATPx = []; 
hiATPy = [];
hiATPstep = [];
[hiATPx, hiATPy, hiATPx_tofit, hiATPy_tofit, hiATPstep, hiATPlabel,...
    hiATPmu, hiATPsigma] = ...
    gatherToCell_Two(hiATPset, stepUpTime, afterstep, ...
    AllRxns_trim, ActTime_trim, PolData_trim, TONORM);

%put all lo-ATP rxns into a single cell array
loATPset = {'sU14'};
loATPx = [];
loATPy = [];
[loATPx, loATPy, loATPx_tofit, loATPy_tofit,loATPstep, loATPlabel, ...
    loATPmu, loATPsigma] = ...
    gatherToCell_Two(loATPset, [nan], afterstep, ... %pass nan instead of stepDownTime b/c fitting only one reaction (=control), which was never stepped-down
    AllRxns_trim, ActTime_trim, PolData_trim, TONORM);

%% fit 

% global fit all high-ATP rxns: sD10, sU1-9 (poststep).
% period same for all rxns, let phase vary
[hiATPfit, hiATPfitdata, hiATPresnorm, hiATPjacobian] = ...
    fitSinusoid_Jac(hiATPx_tofit, hiATPy_tofit, ...
    'period',24,'periodLB',22,'periodUB',26);

% global fit all lo-ATP rxns: sU10, sD1-9 (poststep).
% period same for all rxns, let phase vary
[loATPfit, loATPfitdata, loATPresnorm, loATPjacobian] = ...
    fitSinusoid_Jac(loATPx_tofit, loATPy_tofit, ...
    'period',24,'periodLB',22,'periodUB',26);

%% go through fits to set period the same value for all rxns in a set
for r=1:numel(hiATPfit)
    hiATPfit{r}(2) = hiATPfit{1}(2); %use same period for all fits
end

for r=1:numel(loATPfit)
    loATPfit{r}(2) = loATPfit{1}(2); %use same period for all fits
end

%% turn params into the corresponding 31x1 matrices
% must enclose cell array into another cell array {{[],[],...,[]}}
loATPfit_mat = convCAParToVec({loATPfit});
hiATPfit_mat = convCAParToVec({hiATPfit});

%% test that these new params work in fits
if TOTEST == 1 
    %USE FUNCTION THAT SUPPORTS RESAMPLED DATA
    testParSet_Resampled(hiATPfit_mat,loATPfit_mat,npBootInd);
end


%% compute L and D fns from fit params
  %best way to do this:
    %(1) compute phase of uperturbed sU rxn at step time 
    %     = phase of loATP control at stepUptime
    %     = (2*pi*stepUpTime/loATP_per) - phase_of_control_rxn_in_loATP
    %
    %(2) compute 'perturbed' phase of sU rxn at step time = 
    %     = phase of 'perturbed' sU rxn (in hiATP conditions) at stepUptime
    %     = (2*pi*stepUpTime/hiATP_per) - phase_of_sU_rxn_in_hiATP
    %
    %(3) subtract the two

for n=1:numnew
     
    %hiATP period
    hiATP_fitper = hiATPfit_mat(n,end);
        
    %loATP period
    loATP_fitper = loATPfit_mat(n,end);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Deal with sD rxns %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compute_D_fun = 0;
    if compute_D_fun == 1
        %need phases of hiATP control to compute unperturbed phases 
        % == reference phases
        %did this in two sets, so rxns 1-7 and 8-14 have different
        %reference phases
        hiATP_control_fitphase(1:7) = hiATPfit_mat(n,3); 
        hiATP_control_fitphase(8:14) = hiATPfit_mat(n,6);
        
        %sD phase = phase of hiATP control at stepDownTime
        sDPhase(n,:) = ...
            (2*pi*stepDownTime'/hiATP_fitper) - hiATP_control_fitphase;
        
        %sD phase perturbed = phase of sD rxns at stepDownTime
        %(i.e., use loATP period and phase of fit rxns, not control)
        sDPhase_perturbed(n,:) = ...
            (2*pi*stepDownTime'/loATP_fitper) - loATPfit_mat(n,(1:14)*3);
        
        sDPhaseShift(n,:) = sDPhase - sDPhase_perturbed;
    else
        sDPhase = [];
        sDPhase_perturbed = [];
        sDPhaseShift = [];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Deal with sU rxns %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compute_L_fun = 1;
    if compute_L_fun == 1
        
        %phases of loATP control = references for phase shift calculations
        loATP_control_fitphase = loATPfit_mat(n,end-1); 
        
        %sUPhase = phase of loATP control at stepUpTime
        sUPhase(n,1:13) = ...
            (2*pi*stepUpTime(1:13)'/loATP_fitper) - loATP_control_fitphase;
        
        %sUPhase_perturbed = phase of sU rxn at stepUp time
        sUPhase_perturbed(n,1:13) = ...
            (2*pi*stepUpTime(1:13)'/hiATP_fitper) - hiATPfit_mat(n,(1:13)*3);
        
        %subtract the two
        sUPhaseShift = sUPhase - sUPhase_perturbed;
    else
        sUPhase = [];
        sUPhase_perturbed = [];
        sUPhaseShift = [];
    end
end

%% make L and D periodic
 
%convert to circadian time
z=1; %(2*pi)/(2*pi); %conversion factor

%make sure step phases are in [0, 2*pi]
if max(sUPhase) > 2*pi*0.25
    sUPhase = sUPhase - 2*pi;
    disp('fixed sU ph a');
elseif min(sUPhase) < -2*pi*1.25
    sUPhase = sUPhase + 2*pi;
    disp('fixed sU ph b');
end

if max(sDPhase) > 2*pi*0.25
    sDPhase = sDPhase - 2*pi;
    disp('fixed sD ph a');
elseif min(sDPhase) < -2*pi*1.25
    sDPhase = sDPhase + 2*pi;
    disp('fixed sD ph b');
end

%arrange for double plotting
stepUpPhase = [z*sUPhase z*(sUPhase+2*pi)];
stepDownPhase = [z*sDPhase z*(sDPhase+2*pi)];
stepUpPhaseShift = [z*sUPhaseShift z*sUPhaseShift];
stepDownPhaseShift = [z*sDPhaseShift z*sDPhaseShift];

fTestFns = figure();
plot(stepUpPhase,stepUpPhaseShift,'bs-');
hold on;
%plot(stepDownPhase,stepDownPhaseShift,'rs-');

%% wrap around -- better; think these next 6 lines do a better job than the 
% entire complicated procedure below that compares things to lin fits
% stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) = ...
%     stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) + z*2*pi;
stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) = ...
    stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) + z*2*pi;
if (min(stepUpPhaseShift > 0))
    stepUpPhaseShift = stepUpPhaseShift - 2*pi;
end

%sort step funs after double plotting
%[~,sdOrd] = sort(stepDownPhase);
[~,suOrd] = sort(stepUpPhase);
stepUpPhase = stepUpPhase(suOrd);
stepUpPhaseShift = stepUpPhaseShift(suOrd);
% stepDownPhase = stepDownPhase(sdOrd);
% stepDownPhaseShift = stepDownPhaseShift(sdOrd);

fTestFns2 = figure();
plot(stepUpPhase,stepUpPhaseShift,'bs-');
hold on;
%plot(stepDownPhase,stepDownPhaseShift,'rs-');

%% do linear fits
for n=1:numnew
    xx = z*(-2.5*pi:0.02:2.5*pi);
    
    llo = 4*(2*pi/24);
    lhi = 16*(2*pi/24);
    lbreak = llo; %this is a value in radians, not an actual datapoint
    llo_ix = find(stepUpPhase(n,:) > llo,1,'first');
    lhi_ix = find(stepUpPhase(n,:) < lhi,1,'last');
    lbreak_ix = find(stepUpPhase(n,:) > lbreak,1,'first');
    
    LInd = [llo_ix lhi_ix lbreak_ix];
    
    %find indices for fit -- LInd, DInd aren't for sorted
    sUlinfit(n,:) = fitToLine(stepUpPhase(n,LInd(1):LInd(2)),...
        stepUpPhaseShift(n,LInd(1):LInd(2)));
    sDlinfit = [];
    
    %get [0, 2pi] right, then replicate
    sULineX(n,1:numel(xx)) = xx;
    upBreak(n) = lbreak; 
    
    sULineY(n,xx > upBreak(n) & xx <= upBreak(n)+2*pi*z) = ...
        polyval(sUlinfit(n,[1 2]),...
        xx(xx > upBreak(n) & xx <= upBreak(n)+2*pi*z));
    
    sULineY(n,xx <= upBreak(n) & xx > upBreak(n)-2*pi*z) = ...
        polyval(sUlinfit(n,[1 2])+[0 2*pi*z*sUlinfit(n,1)],...
        xx(xx <= upBreak(n) & xx > upBreak(n)-2*pi*z));
    
    sULineY(n,xx <= upBreak(n)-2*pi & xx > upBreak(n)-4*pi*z) = ...
        polyval(sUlinfit(n,[1 2])+2*[0 2*pi*z*sUlinfit(n,1)],...
        xx(xx <= upBreak(n)-2*pi & xx > upBreak(n)-4*pi*z));
    
    sULineY(n,xx > upBreak(n)+2*pi*z & xx <= upBreak(n)+4*pi*z) = ...
        polyval(sUlinfit(n,[1 2])-[0 2*pi*z*sUlinfit(n,1)],...
        xx(xx > upBreak(n)+2*pi*z & xx <= upBreak(n)+4*pi*z));
    
    sULineY(n, xx > upBreak(n) + 4*pi*z) = ... %just added
        polyval(sUlinfit(n,[1 2])-[0 4*pi*z*sUlinfit(n,1)],...
        xx(xx > upBreak(n) + 4*pi*z));
    
    %get [0, 2pi] right, then replicate
    sDLineX = [];
    downBreak = [];
    sDLineY = [];   
end

%% package outputs (end of required steps)
PARSET_hiATP = hiATPfit_mat;
PARSET_loATP = loATPfit_mat;
T_hiATP = hiATPfit_mat(n,end);
T_loATP = loATPfit_mat(n,end);

STEPOUT = {stepUpPhase, stepUpPhaseShift, T_hiATP, sULineX, sULineY, ...
    sUlinfit, upBreak, ...
    stepDownPhase, stepDownPhaseShift, T_loATP, sDLineX, sDLineY, ...
    sDlinfit, downBreak, npBootInd};

%% plot
if TOPLOT_STEPFUNS == 1

z=(1/(2*pi));
    
fBootStep=figure();
for n=1:numnew
    [~,uind]=sort(stepUpPhase(n,:));
    pUp=plot(z*stepUpPhase(n,uind),z*stepUpPhaseShift(n,uind),'bs-',...
        'markerfacecolor','b','markersize',4,'linewidth',2);
    hold on;
    pUpLin=[];%
    pUpLin=plot(z*sULineX(n,:),z*sULineY(n,:),'b-','linewidth',0.5);
    plot(z*stepUpPhase(n,LInd([1 2])),z*stepUpPhaseShift(n,LInd([1 2])),...
        'ko','markerfacecolor','k','markersize',4,'linewidth',0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %comment this out if you want to plot 
% %    stepUp/stepDown on the same plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(gca,'xlim',z*[-2*pi 2*pi],'xtick',z*[-2*pi:pi:4*pi],...
%     'ylim',z*[-3*pi 3*pi],'ytick',z*[-pi:0.5*pi:pi]);
% xlabel('step time (CT hours)');
% ylabel('phase shift (CT hours)');
% 
% if TOEXPORT_fBootStepUp == 1
%     export_fig([getDate('yyyy-mm-dd') '_stepUp_line_boot_' getDate()...
%         '.pdf'],'-cmyk','-painters','-pdf',fBootStepUp);
% end
% 
% fBootStepDown=figure();
% %endcomment  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for n=1:numnew
%     [~,dind]=sort(stepDownPhase(n,:));
%     pDown=plot(z*stepDownPhase(n,dind),z*stepDownPhaseShift(n,dind),'rs-',...
%         'markerfacecolor','r','markersize',4,'linewidth',2);
%     hold on;
%     pDownLin=[];%
%     pDownLin=plot(z*sDLineX(n,:),z*sDLineY(n,:),'r-','linewidth',0.5);
% end
legend([pUp,pUpLin],...
    'L(\theta)','L_{lin}(\theta)',...
    'location','northwest','orientation','horizontal');
legend boxoff;
set(gca,'xlim',z*[-3*pi 3*pi],'xtick',z*[-2*pi:pi:4*pi],...
    'ylim',z*[-0.5*pi 0.5*pi],'ytick',z*[-2*pi:0.25*pi:2*pi]);
xlabel('step phase \theta (rad/2\pi)');
ylabel('phase shift \Delta\theta (rad/2\pi)');
grid off;
set(fBootStep,'units','inches','position',[0 0 8 4]);

if TOEXPORT_fBootStep == 1
export_fig([getDate('yyyy-mm-dd') '_stepUpDown_line_BESTFIT_noboot' getDate()...
    '.pdf'],'-cmyk','-painters','-pdf',fBootStep);
end
end

end

