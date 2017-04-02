%2017-03-30, EL: do one global fit of all polarization data from a
%stepUp/stepDown experiment. No resampling. 
%
% dependencies: export_fig.m from FileExchange for figure export.
%
function [PARSET_hiATP, PARSET_loATP, STEPOUT] = ...
    bestFitNonParBoot
%Return a set of best fit parameters for a set of resampled loATP and hiATP
%datasets, as well as step functions (and their linearizations) that result
%from this best fit parameter set.

%% load raw data, step times; set normalization/fit parameters
INDIR = ['saved_data'];
INFILE = [INDIR '/' '2017-03-01_stepUpDown_PlateReader.xlsx'];
INSHEETS = {'Sheet3'};

TORESAMPLE = 0; %set this to 0; no resampling, fit all data together

TOTEST = 0; %plot how well the fits fit the data?
TOPLOT_STEPFUNS = 1; %plot step functions?
TOEXPORT_fBootStep = 0; %export figure?

% script to load plate reader data
[AllRxns, AllWells, ActTime, PolData] = loadPlateReader(INFILE,INSHEETS);

%'stepTimes' sheet has stepTimes relative to t=0 on plate reader
INFILE_STEP = [INDIR '/' '2017-02-28_stepUpDown_schedule.xlsx'];
[stepTimes,~,~] = xlsread(INFILE_STEP,'stepTimes','B13:B44');
stepUpTime = stepTimes([1:7 9:15]); 
stepDownTime = [stepTimes(16+[1:7 9:15])];

afterstep = 2; %how much time to wait after step to start fitting
TONORM = 1; %normalize to mean0,std100 or to keep as raw data for fits/plots
TONORMTRIM = 0; %don't normalize while trimming
%%
allRxnTypes = {'sD1','sD2','sD3','sD4','sD5','sD6','sD7','sD8','sD9',...
    'sD10','sD11','sD12','sD13','sD14','sD15','sD16',...
               'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9',...
    'sU10','sU11','sU12','sU13','sU14','sU15','sU16'};...};
         
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
hiATPset = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7',...
            'sU9','sU10','sU11','sU12','sU13','sU14','sU15',...
            'sD8','sD16'}; 
        %sU1-7 have been stepped-up into ATP. 
        %sD8=hiATP MM before ADP addn for step-down (sD1-7).
hiATPx = []; 
hiATPy = [];
hiATPstep = [];
[hiATPx, hiATPy, hiATPx_tofit, hiATPy_tofit, hiATPstep, hiATPlabel,...
    hiATPmu, hiATPsigma] = ...
    gatherToCell(hiATPset, stepUpTime, afterstep, ...
    AllRxns_trim, ActTime_trim, PolData_trim, TONORM);

%put all lo-ATP rxns into a single cell array
loATPset = {'sD1','sD2','sD3','sD4','sD5','sD6','sD7',...
            'sD9','sD10','sD11','sD12','sD13','sD14','sD15',...
            'sU8','sU16'};
loATPx = [];
loATPy = [];
[loATPx, loATPy, loATPx_tofit, loATPy_tofit,loATPstep, loATPlabel, ...
    loATPmu, loATPsigma] = ...
    gatherToCell(loATPset, stepDownTime, afterstep, ...
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

% disp('loATPfit_mat ->');
% disp(loATPfit_mat);
% 
% disp('hiATPfit_mat ->');
% disp(hiATPfit_mat);

%% test that these new params work in fits
if TOTEST == 1 
    %USE FUNCTION THAT SUPPORTS RESAMPLED DATA
    testParSet_Resampled(hiATPfit_mat,loATPfit_mat,npBootInd);
end

%% compute L and D fns from the params
for n=1:numnew
    %compute step phases (use control rxn period (par(n=43)) and phase
    %(par(n,30))
    %recall that sUTime and sDTime have 'nan' as last entry
    
    %compute phase of stepDown by backtracking from phases of control
    %reactions
    sDPhase(n,1:7) = (2*pi*stepDownTime(1:7)/hiATPfit_mat(n,end)) -...
        hiATPfit_mat(n,15*3);
    sDPhase(n,8:14) = (2*pi*stepDownTime(8:14)/hiATPfit_mat(n,end)) -...
        hiATPfit_mat(n,16*3);
    
    %
    sUPhase(n,1:7) = (2*pi*stepUpTime(1:7)/loATPfit_mat(n,end)) -...
        loATPfit_mat(n,15*3);
    sUPhase(n,8:14) = (2*pi*stepUpTime(8:14)/loATPfit_mat(n,end)) -...
        loATPfit_mat(n,16*3);
    
    %compute phase shifts
    sUPhaseShift(n,1:7) = hiATPfit_mat(n,3*(1:7)) - loATPfit_mat(n,3*15);
    sUPhaseShift(n,8:14) = hiATPfit_mat(n,3*(8:14)) - loATPfit_mat(n,3*16);
    
    sDPhaseShift(n,1:7) = loATPfit_mat(n,3*(1:7)) - hiATPfit_mat(n,3*15);
    sDPhaseShift(n,8:14) = loATPfit_mat(n,3*(8:14)) - hiATPfit_mat(n,3*16);
end


%% make L and D periodic
 
%convert to circadian time
z=1; %(2*pi)/(2*pi); %conversion factor

%make sure step phases are in [0, 2*pi]
if max(sUPhase) > 2*pi*0.25
    sUPhase = sUPhase - 2*pi;
    disp('fixed sU ph a');
elseif min(sUPhase) < -2*pi*1.25 || max(sUPhase) < 2*pi
    sUPhase = sUPhase + 2*pi;
    disp('fixed sU ph b');
end

if max(sDPhase) > 2*pi*0.25
    sDPhase = sDPhase - 2*pi;
    disp('fixed sD ph a');
elseif min(sDPhase) < -2*pi*1.25 || max(sDPhase) < 2*pi
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
plot(stepDownPhase,stepDownPhaseShift,'rs-');

%% wrap around -- better; think these next 6 lines do a better job than the 
% entire complicated procedure below that compares things to lin fits
stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) = ...
    stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) + z*2*pi;
% stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) = ...
%     stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) + z*2*pi;
if (min(stepUpPhaseShift > 0))
    stepUpPhaseShift = stepUpPhaseShift - 2*pi;
end

%sort step funs after double plotting
[~,sdOrd] = sort(stepDownPhase);
[~,suOrd] = sort(stepUpPhase);
stepUpPhase = stepUpPhase(suOrd);
stepUpPhaseShift = stepUpPhaseShift(suOrd);
stepDownPhase = stepDownPhase(sdOrd);
stepDownPhaseShift = stepDownPhaseShift(sdOrd);

fTestFns2 = figure();
plot(stepUpPhase,stepUpPhaseShift,'bs-');
hold on;
plot(stepDownPhase,stepDownPhaseShift,'rs-');

%% do linear fits
%
xx=[];
sULineX = [];
sDLineY = [];
upBreak = [];
sUlinfit = [];

sDLineX = [];
sULineY = [];
downBreak = [];
sDlinfit = [];

for n=1:numnew
      DInd = [1 7 1]; %was [2 5 2] %[linLo linHi BreakPt]; breakPt = first pt on new curve
      LInd = [5 14 3]; %include shallower points
      LInd = [9 14 9]; %use only steep region
      
      xx = z*(-2.5*pi:0.1:2.5*pi); %was 0.1
      
%       %find indices for fit -- LInd, DInd aren't for sorted 
%            
      sUlinfit(n,:) = fitToLine(stepUpPhase(n,LInd(1):LInd(2)),...
          stepUpPhaseShift(n,LInd(1):LInd(2)));
      sDlinfit(n,:) = fitToLine(stepDownPhase(n,DInd(1):DInd(2)),...
          stepDownPhaseShift(n,DInd(1):DInd(2)));
%       
      %get [0, 2pi] right, then replicate
      sULineX(n,1:numel(xx)) = xx;
      upBreak(n) = stepUpPhase(n,LInd(3))-0.75; %used to be -0.1
      
      sULineY(n,xx > upBreak(n) & xx <= upBreak(n)+2*pi*z) = ...
          polyval(sUlinfit(n,[1 2]),...
                  xx(xx > upBreak(n) & xx <= upBreak(n)+2*pi*z));
            
      sULineY(n,xx <= upBreak(n) & xx > upBreak(n)-2*pi*z) = ...
          polyval(sUlinfit(n,[1 2])+[0 2*pi*z*sUlinfit(n,1)],...
          xx(xx <= upBreak(n) & xx > upBreak(n)-2*pi*z));
      
      sULineY(n,xx > upBreak(n)+2*pi*z & xx <= upBreak(n)+4*pi*z) = ... 
              polyval(sUlinfit(n,[1 2])-[0 2*pi*z*sUlinfit(n,1)],...
                  xx(xx > upBreak(n)+2*pi*z & xx <= upBreak(n)+4*pi*z));
      
      sULineY(n, xx > upBreak(n) + 4*pi*z) = ... %just added
              polyval(sUlinfit(n,[1 2])-[0 4*pi*z*sUlinfit(n,1)],...
                  xx(xx > upBreak(n) + 4*pi*z));
    
        
      %get [0, 2pi] right, then replicate
      sDLineX(n,1:numel(xx)) = xx;
      downBreak(n) = stepDownPhase(n,DInd(3))-0.1; %kluge to get breakpt just after step
              
      sDLineY(n,xx > downBreak(n) & xx <= downBreak(n)+2*pi*z) = ...
          polyval(sDlinfit(n,[1 2]),...
                  xx(xx > downBreak(n) & xx <= downBreak(n)+2*pi*z));
              
      sDLineY(n,xx <= downBreak(n)) = ...
          polyval(sDlinfit(n,[1 2])+[0 2*pi*z*sDlinfit(n,1)],...
          xx(xx <= downBreak(n)));
      
      sDLineY(n,xx > downBreak(n)+2*pi*z & xx <= downBreak(n)+4*pi*z) = ... 
              polyval(sDlinfit(n,[1 2])-[0 2*pi*z*sDlinfit(n,1)],...
                  xx(xx > downBreak(n)+2*pi*z & xx <= downBreak(n)+4*pi*z));
              
      sDLineY(n, xx > downBreak(n)+4*pi*z) = ...
              polyval(sDlinfit(n,[1 2]) - 2*[0 2*pi*z*sDlinfit(n,1)],...
                  xx(xx > downBreak(n) + 4*pi*z));
                    
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
TOPLOT_STEPFUNS = 1;
if TOPLOT_STEPFUNS == 1

z=(1/(2*pi));
    
fBootStep=figure();
for n=1:numnew
    [~,uind]=sort(stepUpPhase(n,:));
    pUp=plot(z*stepUpPhase(n,uind),z*stepUpPhaseShift(n,uind),'bs-',...
        'markerfacecolor','b','markersize',4,'linewidth',2);
    hold on;
    plot(z*stepUpPhase(n,LInd([1 2])),z*stepUpPhaseShift(n,LInd([1 2])),...
        'yo','markerfacecolor','y','markersize',4,'linewidth',0.5);
   
    pUpLin=[];%
    pUpLin=plot(z*sULineX(n,:),z*sULineY(n,:),'b-','linewidth',0.5);
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

for n=1:numnew
    [~,dind]=sort(stepDownPhase(n,:));
    pDown=plot(z*stepDownPhase(n,dind),z*stepDownPhaseShift(n,dind),'rs-',...
        'markerfacecolor','r','markersize',4,'linewidth',2);
    hold on;
    plot(z*stepDownPhase(n,DInd([1 2])),z*stepDownPhaseShift(n,DInd([1 2])),...
        'ko','markerfacecolor','k','markersize',4,'linewidth',0.5);
    pDownLin=[];%
    pDownLin=plot(z*sDLineX(n,:),z*sDLineY(n,:),'r-','linewidth',0.5);
end
legend([pDown,pDownLin,pUp,pUpLin],...
    'D(\theta)','D_{lin}(\theta)',...
    'L(\theta)','L_{lin}(\theta)',...
    'location','northwest','orientation','horizontal');
legend boxoff;
set(gca,'xlim',z*[-2*pi 2*pi],'xtick',z*[-2*pi:pi:4*pi],...
    'ylim',z*[-pi pi],'ytick',z*[-2*pi:0.25*pi:2*pi]);
xlabel('step phase \theta (rad/2\pi)');
ylabel('phase shift \Delta\theta (rad/2\pi)');
grid off;
set(fBootStep,'units','inches','position',[0 0 8 4]);

if TOEXPORT_fBootStep == 1
export_fig([getDate('yyyy-mm-dd') '_stepUpDown_line_BESTFIT_noboot_w_shallower_' getDate()...
    '.pdf'],'-cmyk','-painters','-pdf',fBootStep);
end
end

end

