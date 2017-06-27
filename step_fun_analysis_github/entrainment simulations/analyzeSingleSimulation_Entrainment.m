%2017-03-19, EL: used to make Fig. 4B, Fig. 4-figSup1(A)
%run entrainment simulations for multiple-light dark cycles starting from
%different initial conditions.
%
% This script runs simulations by using step-up and step-down functions in
% the following file: STEPFUN_BOOT=[FILENAME].
%
% dependencies: export_fig.m from FileExchange to export figures.

clear all;
close all;
clc;

INDIR=['.'];
cd(INDIR);

%export figs?
toexp_fDiffTau = 0; %(1=yes) simulations for different tau values from same starting phase
toexp_fDiffPhases = 0; %(1=yes) simulations for same tau value for different starting phases

%% set simulation parameters here
TDRIVE=24.0; %what driving period?
dutyfrac=[8 10 12 14]/24; %what day lengths?
numcyc = 1:6; % how many LD cycles to simulate?
TEND = (max(numcyc)+0.1)*TDRIVE; % how long to run simulation?

%what is the initial phase for simulations?
startphase = [0 0.2 0.4 0.8];

%which step function number to use out of bootstrapped ones in the file?
%must provide one number < total number step functions in file!
boot = [1];
assert((numel(boot)) == 1); %analyze one a at a time

%% load step function
%NOTE! must be in 2*pi units
STEPFUN_BOOT = ['../helper functions and shared files/'...
    '2017-06-05_widefit_mergedStepFuns_11.30.54.mat'];
load([STEPFUN_BOOT]);

%% run simulations with different L and D funs
for b=boot
       
    TLIGHT=T_hiATP_mix(b);
    TDARK=T_loATP_mix(b);
    STEPUP.phase = up_mix{b}.phase;
    STEPUP.phaseShift = up_mix{b}.phaseShift;
    STEPDOWN.phase = down_mix{b}.phase;
    STEPDOWN.phaseShift = down_mix{b}.phaseShift;

    STEPUP.linPhase = up_mix{b}.linPhase;
    STEPUP.linPhaseShift = up_mix{b}.linPhaseShift;
    STEPDOWN.linPhase = down_mix{b}.linPhase;
    STEPDOWN.linPhaseShift = down_mix{b}.linPhaseShift;
   
    %time this
    tictime=tic;
    disp(['Starting boot ' num2str(b) '...']);
    for df=1:numel(dutyfrac)
        for s=1:numel(startphase)
           
            try
                %pass step funs to be interpolayed to this function
                [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
                    drivePhaseOscilStepFn_WithBoot_Fast_Two(TDRIVE, TLIGHT, TDARK,...
                    STEPUP, STEPDOWN, dutyfrac(df), TEND, ...
                    numcyc, startphase(s));
                
                peakT(df, 1:numel(numcyc),s) = (0.5-(DAWNPHASE+DAWNPHASESHIFT)')./(1/TLIGHT);
                dawnphase(df,1:numel(numcyc),s) = DAWNPHASE;
                duskphase(df,1:numel(numcyc),s) = DUSKPHASE;
                duskphaseshift(df,1:numel(numcyc),s) = DUSKPHASESHIFT;
                dawnphaseshift(df,1:numel(numcyc),s) = DAWNPHASESHIFT;
            catch err
                disp(['boot ' num2str(b)  ...
                    ', theta_0 ' num2str(startphase(s)) ...
                    ', df ' num2str(dutyfrac(df)) ': simulation failed']);
                peakT(df, 1:numel(numcyc),s) = nan;
                dawnphase(df,1:numel(numcyc),s) = nan;
                duskphase(df,1:numel(numcyc),s) = nan;
                duskphaseshift(df,1:numel(numcyc),s) = nan;
                dawnphaseshift(df,1:numel(numcyc),s) = nan;
            end
        end
    end
    
    elapsedtime=toc(tictime)/60;
    disp(['Boot no. ' num2str(b) ...
        ' completed in ' num2str(elapsedtime) ' min.']);
end

%% plot phase as a fn of num. LD cycles for multiple starting points
fDiffPhases = figure();

plotphase = dawnphase+dawnphaseshift;
plotphase = mod(plotphase,1);
plotphase(plotphase > 0.8) = plotphase(plotphase > 0.8) - 1;

p=[];
s=[];
lab=[];

for df=1:numel(dutyfrac)
    figPhases(df)=figure;
    colors = hsv(numel(startphase));
    p=[];
    for s=1:numel(startphase)
        p(s) = plot([0 numcyc],[startphase(s) plotphase(df,1:numel(numcyc),s)],...
            's','markersize',4,'markerfacecolor',colors(s,:),...
            'markeredgecolor','none');
        hold on;
        plot([0 numcyc],[startphase(s) plotphase(df,1:numel(numcyc),s)],...
            '-','color',colors(s,:),'linewidth',1);
        lab{s} = ['\theta_0=' num2str(startphase(s),'%2.1f')];
    end
    if df>0
        legend(p,lab,'location','northeast');
    end
    legend boxoff;
    xlabel('days');
    ylabel('oscillator phase at dawn (rad/2\pi)');
    title(['LD' num2str(TDRIVE*dutyfrac(df)) ':' ...
        num2str(TDRIVE*(1-dutyfrac(df)))]);
    set(gca,'xlim',[0 max(numcyc)],'ylim',[-0.2 0.8], ...
    'xtick',0:1:max(numcyc),'ytick',[-0.2:0.2:1],'fontsize',12);
    grid off;
    set(figPhases(df),'units','inches','position',[0 0 3 3]);

if toexp_fDiffPhases == 1
    export_fig([pwd '/' getDate('yyyy-mm-dd') '_PO_diff-theta0_tau-'... ...
        num2str(dutyfrac(df)*24) '_' getDate('HH.MM.SS')],...
        '-cmyk','-painters','-pdf',figPhases(df));
end

end

%% plot phase as a fn of num. LD cycles for multiple tau's
fDiffTau = figure();

plotphase = dawnphase+dawnphaseshift;
plotphase = mod(plotphase,1);
plotphase(plotphase > 0.8) = plotphase(plotphase > 0.8) - 1;

p=[];
s=[];
lab=[];

for s=1:numel(startphase)
    figTau(s)=figure;
    colors = hsv(numel(dutyfrac));
    p=[];
    lab=[];
    for df=1:numel(dutyfrac)
        p(df) = plot([0 numcyc],[startphase(s) plotphase(df,1:numel(numcyc),s)],...
            's','markersize',4,'markerfacecolor',colors(df,:),...
            'markeredgecolor','none');
        hold on;
        plot([0 numcyc],[startphase(s) plotphase(df,1:numel(numcyc),s)],...
            '-','color',colors(df,:),'linewidth',1);
        lab{df} = ['\tau = ' num2str(TDRIVE*dutyfrac(df)) ' h'];
    end
    if s>0
        legend(p,lab,'location','northeast');
    end
    legend boxoff;
    xlabel('days');
    ylabel('oscillator phase at dawn (rad/2\pi)');
    %title(['\theta_0=' num2str(startphase(s))]);
    set(gca,'xlim',[0 max(numcyc)],'ylim',[-0.2 0.8], ...
    'xtick',0:1:max(numcyc),'ytick',[-0.2:0.2:1],'fontsize',12);
    grid off;
    
    set(figTau(s),'units','inches','position',[0 0 3 3]);

if toexp_fDiffTau == 1
    export_fig([pwd '/' getDate('yyyy-mm-dd') '_PO_diff-taus_theta0-'...
        num2str(startphase(s)) '_' getDate('HH.MM.SS')],...
        '-cmyk','-painters','-pdf',figTau(s));
end

end