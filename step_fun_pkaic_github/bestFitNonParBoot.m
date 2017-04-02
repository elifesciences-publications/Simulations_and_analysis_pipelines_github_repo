%2017-03-29, EL: perform one non-parametric bootstrap resampling of
%stepup/down data and do a global fit, as described in Computational
%Methods. 
%
% Change function parameters below to determine whether you want to turn
% resampling on/off, whether to plot the best fit to resampled data and 
% whether to plot the step functions at the end of the run.
%
% This function can be run as a function or one Matlab cell at a time (to
% store variables into workspace).
%
% dependencies: export_fig.m from FileExchange for figure export. 

function [PARSET_hiATP, PARSET_loATP, STEPOUT] = ...
    bestFitNonParBoot
%Return a set of best fit parameters for a set of resampled loATP and hiATP
%datasets, as well as step functions (and their linearizations) that result
%from this best fit parameter set.

%% function parameters
TOTEST = 0; %to run function that plots resampled data (one bootstrap) + best fit to it?
TORESAMPLE = 1; %if 1, sample all data with replacement ; if 0, don't bootstrap -- use entire dataset
TOPLOT_STEPFUNS = 0; %plot step functions for this sample of data?
TOEXPORT_fBootStep = 0; %export step fun figure?

%% load raw data, step times; set normalization/fit parameters
INDIR = ['saved_data'];
INFILE = [INDIR '/2016-08-25_all_quants.xlsx'];

% script to load SDS-PAGE quants of KaiC phosphorylation data
[AllRxns, TubeTime, ActTime, pSpT, pT, pS, pU, pNet] = loadQuants(INFILE); 
stepUpTime = [0 2 4 6 8 12 16 20 24 nan];
stepDownTime = [0 4 8 12 14 16 18 20 24 nan];

afterstep = 16; %how much time to wait after step to start fitting
TONORM = 1; %normalize to mean0,std100 or to keep as raw data for fits/plots
TONORMTRIM = 0; %don't normalize while trimming

allRxnTypes = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9','sU10',...
               'sD1','sD2','sD3','sD4','sD5','sD6','sD7','sD8','sD9','sD10'};
         
[ActTime_trim, pNet_trim, AllRxns_trim] = ...
    gatherToMat(allRxnTypes, [stepUpTime stepDownTime], ...
    afterstep, AllRxns, ActTime, pNet, TONORMTRIM); %do not normalize here
AllRxns_trim = AllRxns_trim'; %must transpose

%% draw a random set of datapoints with replacement
numpts = numel(pNet_trim); %total num. pts

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
pNet_trim = pNet_trim(npBootInd);
ActTime_trim = ActTime_trim(npBootInd);

%% gather new datapoints into cell arrays
%put all high-ATP rxns into a single cell array
hiATPset = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9','sD10'};
hiATPx = []; 
hiATPy = [];
hiATPstep = [];
[hiATPx, hiATPy, hiATPx_tofit, hiATPy_tofit, hiATPstep, hiATPlabel,...
    hiATPmu, hiATPsigma] = ...
    gatherToCell(hiATPset, stepUpTime, afterstep, ...
    AllRxns_trim, ActTime_trim, pNet_trim, TONORM);

%put all lo-ATP rxns into a single cell array
loATPset = {'sD1','sD2','sD3','sD4','sD5','sD6','sD7','sD8','sD9','sU10'};
loATPx = [];
loATPy = [];
[loATPx, loATPy, loATPx_tofit, loATPy_tofit,loATPstep, loATPlabel, ...
    loATPmu, loATPsigma] = ...
    gatherToCell(loATPset, stepDownTime, afterstep, ...
    AllRxns_trim, ActTime_trim, pNet_trim, TONORM);

%% fit 

% global fit all high-ATP rxns: sD10, sU1-9 (poststep).
% period same for all rxns, let phase vary
[hiATPfit, hiATPfitdata, hiATPresnorm, hiATPjacobian] = ...
    fitSinusoid_Jac(hiATPx_tofit, hiATPy_tofit, 'period',24,'periodLB',22,'periodUB',26);

% global fit all lo-ATP rxns: sU10, sD1-9 (poststep).
% period same for all rxns, let phase vary
[loATPfit, loATPfitdata, loATPresnorm, loATPjacobian] = ...
    fitSinusoid_Jac(loATPx_tofit, loATPy_tofit, 'period',24,'periodLB',22,'periodUB',26);

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

%return;
%% compute L and D fns from the params
for n=1:numnew
    %compute step phases (use control rxn period (par(n=31)) and phase
    %(par(n,30))
    %recall that sUTime and sDTime have 'nan' as last entry
    sUPhase(n,:) = (2*pi*stepUpTime(1:9)/hiATPfit_mat(n,end)) - hiATPfit_mat(n,30);
    sDPhase(n,:) = (2*pi*stepDownTime(1:9)/loATPfit_mat(n,end)) - loATPfit_mat(n,30);
    
    %compute phase shifts
    sUPhaseShift(n,:) = hiATPfit_mat(n,3:3:27) - hiATPfit_mat(n,30);
    sDPhaseShift(n,:) = loATPfit_mat(n,3:3:27) - loATPfit_mat(n,30);
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

%wrap around
% stepDownPhaseShift(stepDownPhaseShift < -(z*pi*0.5)) = ...
%     stepDownPhaseShift(stepDownPhaseShift < -(z*pi*0.5)) + z*2*pi;
% stepDownPhaseShift(stepDownPhaseShift > z*1.2*pi) = ...
%     stepDownPhaseShift(stepDownPhaseShift > z*1.2*pi) - z*2*pi;
% stepUpPhaseShift(stepUpPhaseShift < -(z*pi*1.1)) = ...
%     stepUpPhaseShift(stepUpPhaseShift < -(z*pi*1.1))+z*2*pi;

%% wrap around -- better; think these next 6 lines do a better job than the 
% entire complicated procedure below that compares things to lin fits
stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) = ...
    stepDownPhaseShift(stepDownPhaseShift < stepDownPhaseShift(1)-0.1) + z*2*pi;
stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) = ...
    stepUpPhaseShift(stepUpPhaseShift < stepUpPhaseShift(6)-0.1) + z*2*pi;
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

%% do linear fits
for n=1:numnew

      DInd = [2 5 1]; %was [2 5 2] %[linLo linHi BreakPt]; breakPt = first pt on new curve
      LInd = [6 10 6];
      
      xx = z*(-2.5*pi:0.02:2.5*pi); %was 0.1
      
      %find indices for fit -- LInd, DInd aren't for sorted 
           
      sUlinfit(n,:) = fitToLine(stepUpPhase(n,LInd(1):LInd(2)),...
          stepUpPhaseShift(n,LInd(1):LInd(2)));
      sDlinfit(n,:) = fitToLine(stepDownPhase(n,DInd(1):DInd(2)),...
          stepDownPhaseShift(n,DInd(1):DInd(2)));
      
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

%% adjust kinks
stepUpPhaseShift(n,:) = adjKinksToLinear(stepUpPhase(n,:),stepUpPhaseShift(n,:),...
                                            sULineX(n,:),sULineY(n,:));
stepDownPhaseShift(n,:) = adjKinksToLinear(stepDownPhase(n,:), stepDownPhaseShift(n,:),...
                                            sDLineX(n,:), sDLineY(n,:));
                                        
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
    pDownLin=plot(z*sDLineX(n,:),z*sDLineY(n,:),'r-','linewidth',0.5);
end
legend([pUp,pDown,pUpLin,pDownLin],...
    'L(\theta)','D(\theta)','L_{lin}(\theta)','D_{lin}(\theta)',...
    'location','northwest','orientation','horizontal');
legend boxoff;
set(gca,'xlim',z*[-2*pi 2*pi],'xtick',z*[-2*pi:pi:4*pi],...
    'ylim',z*[-2*pi 2*pi],'ytick',z*[-2*pi:1*pi:2*pi]);
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

