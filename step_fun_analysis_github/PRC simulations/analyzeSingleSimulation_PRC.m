%2017-03-28, EL: generates Fig. 5-figSup1
%simulate a phase-resetting curve in the phase oscillator framework
%dependencies: export_fig.m for figure export

clear all;
close all;
clc;

INDIR=['.'];
cd(INDIR);

%% set PRC simulation parameters
TODISP = 0; %display outputs?
LO_PP = 0; % which DPs to mark (LO_PP = phase of earliest DP, etc.);
HI_PP = 24; % HI_PP = last time to give DP
TOEXP_DOUBLEFIG = 1; %save fig?

%%
DPDUR = 12;
TEND = 100;

dutyfrac=[LO_PP:1:HI_PP]/24;
startphase=zeros(size(dutyfrac));

%% set up colors just for pts in linear region
colorInd = dutyfrac >= LO_PP/24 & dutyfrac <= HI_PP/24;
numcolors = sum(colorInd);
colors = zeros(numel(dutyfrac),3);
colors(colorInd,:) = hot(numcolors);
colors = flipud(colors); %so that later DP = darker colors

%% load boostrapped step functions
%NOTE! must be in radians. don't normalize by 2*pi
STEPFUN_BOOT = '../helper functions and shared files/2017-06-05_widefit_mergedStepFuns_11.30.54.mat';
load([STEPFUN_BOOT]);
boot = [2]; %which of the bootstraps or combinations of L and D do you want to run?

%% run PRC simulations
for b=boot
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    STEPUP.phase = up_mix{b}.phase;
    STEPUP.phaseShift = up_mix{b}.phaseShift;
    STEPDOWN.phase = down_mix{b}.phase;
    STEPDOWN.phaseShift = down_mix{b}.phaseShift;
    
    %USE linearized versions
    STEPUP.phase = up_mix{b}.linPhase;
    STEPUP.phaseShift = up_mix{b}.linPhaseShift;
    STEPDOWN.phase = down_mix{b}.linPhase;
    STEPDOWN.phaseShift = down_mix{b}.linPhaseShift;
   
    %time this
    tictime=tic;
    disp(['Starting boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
          
        %catch errors in integration due to kinks
        try
            %pass step funs to be interpolayed to this function
            %DPDUR, DPTIME should both be passed in hrs
            [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                 drivePhaseOscilStepFn_PRC(TEND, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
            dutyfrac(df)*TLIGHT, DPDUR, startphase(df));
            
            peakT(df, b) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
            dawnphase(df,b) = DAWNPHASE;
            duskphase(df,b) = DUSKPHASE;
            duskphaseshift(df,b) = DUSKPHASESHIFT;
            dawnphaseshift(df,b) = DAWNPHASESHIFT;
            pshift(df,b) = (DAWNPHASE+DAWNPHASESHIFT)-...
                (dutyfrac(df)*TLIGHT+DPDUR)./(TLIGHT); %subtract unperturbed
        catch err
            disp(['boot ' num2str(b) ', df ' num2str(dutyfrac(df)) ' simulation failed']);
            peakT(df, b) = nan;
            dawnphase(df,b) = nan;
            duskphase(df,b) = nan;
            duskphaseshift(df,b) = nan;
            dawnphaseshift(df,b) = nan;
            pshift(df,b) = nan;
        end
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end


%% collect variables for plotting
plotphase = pshift;
plotphase = mod(plotphase,1);
plotphase(plotphase > 0.25) = plotphase(plotphase > 0.25) - 1;

%% plot PRC 
fSilico=figure();

subp1=subplot(1,3,1);
for b=boot
    pSilicoA(b)=plot(24*dutyfrac, plotphase(:,b), 'k-',...
        'linewidth',2);
    hold on;
end

%plot with different colors for different taus
for b=boot
    for df=1:numel(dutyfrac)
        if dutyfrac(df) > LO_PP/24 & dutyfrac(df) < HI_PP/24
            pSilico(b)=plot(24*dutyfrac(df), plotphase(df,b), 'o',...
                'markerfacecolor',colors(df,:),'markeredgecolor','k',...
                'markersize',6,'linewidth',0.1);
            hold on;
        end
    end
    %compute slope of simulated curve
    dfRange = dutyfrac > LO_PP/24 & dutyfrac < HI_PP/24;
    simSlope(b,:) = fitToLine(24*dutyfrac(dfRange),plotphase(dfRange,b)');
    disp(['boot ' num2str(b) ',' 10 'sim slope = ' num2str(simSlope(b,1),'%2.3f')]);
    %title(['boostrap no. ' num2str(b)]);
end
colormap(colors);
cbar=colorbar('ytick',[find(colorInd,1,'first'):4:find(colorInd,1,'last')],...
        'yticklabel',24*dutyfrac([find(colorInd,1,'first'):4:find(colorInd,1,'last')]),...
        'ylim',[find(colorInd,1,'first')-0.02 find(colorInd,1,'last')+0.02]);
title(cbar,'t_{DP}');
    
legend boxoff;
xlabel(['dark pulse time ' 't_{DP} (CT hours)']);
ylabel('phase shift (rad/2\pi)');
set(gca,'xlim',[0 24],'ylim',[-0.8 0.4],'xtick',0:6:24,'ytick',[-0.6:0.2:0.6]);
grid off;
%set(fInVitroSilico,'units','inches','position',[0 0 8 8]);

%% plot bootstrapped step functions
%fStepFuns=figure();
subp2=subplot(1,3,[2 3]);
for b=boot
    
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    STEPUP.phase = up_mix{b}.phase;
    STEPUP.phaseShift = up_mix{b}.phaseShift;
    STEPDOWN.phase = down_mix{b}.phase;
    STEPDOWN.phaseShift = down_mix{b}.phaseShift;
    
    %linearized versions
    STEPUP.linPhase = up_mix{b}.linPhase;
    STEPUP.linPhaseShift = up_mix{b}.linPhaseShift;
    STEPUP.linFit = up_mix{b}.linFit;
    STEPDOWN.linPhase = down_mix{b}.linPhase;
    STEPDOWN.linPhaseShift = down_mix{b}.linPhaseShift;
    STEPDOWN.linFit = down_mix{b}.linFit;
    
    z=(1/(2*pi));
    
    figure(fSilico);
    subplot(subp2);
    %double plot interpolated step funs
    LW=0.5;
    UP_LIN.phase = STEPUP.linPhase;
    UP_LIN.phaseShift = STEPUP.linPhaseShift;
    DOWN_LIN.phase = STEPDOWN.linPhase;
    DOWN_LIN.phaseShift = STEPDOWN.linPhaseShift;
    LW=2;
    sUplt=plot([0:0.01:1],-getInVitroPhaseShiftBoot([0:0.01:1],...
        UP_LIN),'b','linewidth',LW);
    hold on;
    plot([0:0.01:1]-1,-getInVitroPhaseShiftBoot([0:0.01:1],...
        UP_LIN),'b','linewidth',LW);
    sDplt=plot([0:0.01:1],-getInVitroPhaseShiftBoot([0:0.01:1],...
        DOWN_LIN),'r','linewidth',LW);
    plot([0:0.01:1]-1,-getInVitroPhaseShiftBoot([0:0.01:1],...
        DOWN_LIN),'r','linewidth',LW);
    
    %mark points used in entrained regime
    linDx = []; linLx = [];
    linDy = []; linLy = [];
    cnt = 1;
    for df=1:numel(dutyfrac)
        if dutyfrac(df) > LO_PP/24 & dutyfrac(df) < HI_PP/24
            plot(wrapVecAround(duskphase(df,b),0.25,1,'gt'),...
                -duskphaseshift(df,b),'o',...
                'markerfacecolor',colors(df,:),'markeredgecolor','k',...
                'linewidth',0.1,'markersize',6);
            plot(wrapVecAround(dawnphase(df,b),0.5,1,'gt'),...
                -dawnphaseshift(df,b),'o',...
                'markerfacecolor',colors(df,:),'markeredgecolor','k',...
                'linewidth',0.1,'markersize',6);
           
            cnt = cnt+1;
        end
    end
    
   
    xlabel('phase (rad/2\pi)');
    ylabel('phase shift (rad/2\pi)');
    set(gca,'xtick',[-1:0.2:1],'ytick',[-0.4:0.1:0.4],...
            'xlim',[-1 1], 'ylim', [-0.4 0.4]);    
    grid off;
    legend([sUplt, sDplt],'step up','step down','location','south');
    legend boxoff;
    
    cbar2=colorbar('ytick',[find(colorInd,1,'first'):4:find(colorInd,1,'last')],...
        'yticklabel',24*dutyfrac([find(colorInd,1,'first'):4:find(colorInd,1,'last')]),...
        'ylim',[find(colorInd,1,'first')-0.02 find(colorInd,1,'last')+0.02]);
    title(cbar2,'t_{DP}');
    
end

set(gcf,'units','inches','position',[0 0 12 4]);
if TOEXP_DOUBLEFIG == 1
    export_fig([getDate('yyyy-mm-dd') '_MixAndMatchNo' num2str(boot) ...
        '_LinBoot_Steps_' getDate() '.pdf'],...
        '-cmyk','-painters',gcf);
end

