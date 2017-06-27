%2017-02-27, EL: modify for merged dataset from Feb5 and Feb23 step down
%measurements
%2017-02-14, EL: modify for plate reader data
%2016-10-31, EL: analyze saved simulations. make plots for manuscript.
%2016-08-22, EL: modify to make plots in phase coordinates for manuscript
%2016-08-16, EL
clear all;
close all;
clc;
cd(['.']);
%%
toexp_VitroSilico = 1; %export figs?
TOSAVE_DATA = 0;
TOLOAD_OLD_SIMULATION_ONE = 1;
TOLOAD_OLD_SIMULATION_TWO = 0;

%% load in vitro data

%compute startphase based on entrained conditions observed in vitro
PHASEFILE='2017-02-15_inVitroPP_Plus1hr_2017-02-15_12.45.23.mat';
PHASEDIR='/Users/E/Documents/Advisers/Rust/Photoperiod/2017-02-05_stepUp_stepDown/analysis';
load([PHASEDIR '/' PHASEFILE]);
%inVitroDawnPhase = (-InVitroFitPhase+0.5*pi)/(2*pi); 
inVitroDawnPhase = inVitroDawnPhase / (2*pi);
inVitroDawnPhase(inVitroDawnPhase > 0.3) = ...
    inVitroDawnPhase(inVitroDawnPhase > 0.3) -1;
inVitroErr = (inVitroPhaseErr)/(2*pi);

%% input simulation parameters -- in the future, these will be stored with the simulations themselves
TDRIVE=24.0;
dutyfrac=[6:1:18]/24;
numcyc = 1:40;
TEND = max(numcyc)*TDRIVE;

%starting phases for each dutyfrac
startphase = interp1(inVitroDuty(4:end), inVitroDawnPhase(4:end), ...
    dutyfrac,'linear','extrap');
startphase = zeros(size(dutyfrac)); %[0];

%how many bootstraps to evaluate?
boot = [1];

%load bootstrapped step funs
% NOTE! must be in 2*pi units
INDIR=['saved_data/'];
INFILE='2017-03-30_stepFun_repeat2_Mar30_test1_2017-03-30_09.53.38.mat';
load([INDIR '/' INFILE]);

%% load results of first simulation here
if TOLOAD_OLD_SIMULATION_ONE == 1
    SIMULATION_BOOT = load('saved_data/2017-03-30_09.56.27_NonParamBoot_I22_HiBrPt_1.mat'); 
    SIM = SIMULATION_BOOT.VAR;
    dawnphase = SIM.dawnphase;
    duskphase = SIM.duskphase;
    duskphaseshift = SIM.duskphaseshift;
    dawnphaseshift = SIM.dawnphaseshift;
    peakT = SIM.peakT;
end

%% collect variables for plotting
%plotphase(df,1:numel(numcyc),s)
plotphase = dawnphase+dawnphaseshift;
plotphase = mod(plotphase,1);
plotphase(plotphase > 0.3) = plotphase(plotphase > 0.3) - 1;

%% select only those tau's where entrained stably
STDEV_LIM = 0.1;

NUMENTR=4; %for nonlinear case
for b=boot
    for df=1:numel(dutyfrac)
        plotphase_stdev(df,b) = std(squeeze(plotphase(df,(end-NUMENTR):end,b)));
        [~,plotphase_stdev(df,b),~] = ...
            circleMean(squeeze(plotphase(df,(end-NUMENTR):end,b)),1,[]);
        if plotphase_stdev(df,b) > STDEV_LIM
            entrained(df,b) = 0;
        else
            entrained(df,b) = 1;
        end
    end
end

%display some entrainment stats
disp(['Total no. sims.: ' num2str(numel(entrained))]);
disp(['No. entrained sims.: ' num2str(sum(sum(entrained)))]);
disp(['Frac entrained: ' num2str(...
    sum(sum(entrained))/...
    numel(entrained)...
    )]);

%specifically for tau between 8 and 16
gf = dutyfrac <= 16/24 & dutyfrac >= 8/24;
disp(['Total no. sims. tau=6-18: ' num2str(numel(entrained(gf,:)))]);
disp(['No. entrained sims. tau=6-18: ' num2str(sum(sum(entrained(gf,:))))]);
disp(['Frac entrained tau=6-18: ' num2str(...
    sum(sum(entrained(gf,:)))/...
    numel(entrained(gf,:))...
    )]);


%% figure out how quickly you entrain
pstdev=[];
numToEnt=[];
NME=3;

for b=boot
    for df=1:numel(dutyfrac)
        for n=1:numel(numcyc)-NME
            %pstdev(df,b,n) = std(squeeze(plotphase(df,n:n+NME,b)));
            [~,pstdev(df,b,n),~] = ...
                circleMean(squeeze(plotphase(df,n:n+NME,b)),1,[]);
        end
        
        %find index of
        nte = find(pstdev(df,b,:) <= STDEV_LIM,1,'first');
        if ~isempty(nte)
            numToEnt(b,df) = numcyc(nte);
        else
            numToEnt(b,df) = nan;
        end
    end
end

%%
%sum all those entrained within 3 cycles, and how many unentr
%for daylengths 6-18 hrs
gf = dutyfrac <= 18/24 & dutyfrac >= 6/24;
gNumToEnt = numToEnt(:,gf);

num1to3 = sum(sum(gNumToEnt == 1 | gNumToEnt == 2 | gNumToEnt == 3));
numNan = sum(sum(isnan(gNumToEnt)));
numNet = numel(gNumToEnt);
fracQuick = num1to3/numNet;
fracQuick_notNan = num1to3/(numNet-numNan);
disp([num2str(fracQuick) ' entrained within 3 cycles']);
disp([num2str(fracQuick_notNan) ' entrained w/in 3cyc out of nonNan']);
% %% check unentrained ones
% fUnentrained=figure();
% for b=1:10%boot
%     for df=1:numel(dutyfrac)
%         if entrained(df,b) == 0
%             plot(numcyc,plotphase(df,:,b));
%             hold all;
%         end
%     end
% end

%% plot entrained phase in silico vs in vitro
fInVitroSilico=figure();

for b=boot
    if sum(entrained(:,b) > 0) %& goodm(b) == 1
        entr = logical(entrained(:,b));
        pSilico(b)=plot(24*dutyfrac(entr), plotphase(entr,end,b),...
            'bs-','markersize',4,...
            'markerfacecolor','b',...
            'linewidth',2);
        hold on;
    end
end

pVitroErr=errorbar(24*inVitroDuty(1:5), inVitroDawnPhase(1:5), ...
    inVitroErr(1:5), 'rs-',...
    'markersize',4,'markerfacecolor','none',...
    'markeredgecolor','none','linewidth',2);
hold on;

legend([pVitroErr,pSilico(1)],'in vitro','simulation','location','south');
legend boxoff;
xlabel('day length \tau (hours)');
ylabel('oscillator phase at dawn (rad/2\pi)');
set(gca,'xlim',[0 24],...
    'ylim',[-1 0],...
    'xtick',0:6:24,'ytick',[-1:0.2:1]);
grid off;
set(fInVitroSilico,'units','inches','position',[0 0 3.75 3]);

if toexp_VitroSilico == 1
    export_fig([pwd '/' getDate('yyyy-mm-dd') '_fInVitroInSilico_'...
        'NonParamBoot' num2str(max(boot)) '_' getDate('HH.MM.SS')],...
        '-cmyk','-painters','-pdf',fInVitroSilico);
end

%% plot avg phase+/-std for entrained tau's for entire boostrapped set
%go through dfs, collect [tau, taustd]
for df=1:numel(dutyfrac)
    tau(df) = dutyfrac(df);
    entr = logical(entrained(df,:));
    goodph = plotphase(df,end,entr);
    goodph = goodph(goodph < 0.33);
    %     taumean(df) = nanmean(goodph);
    %     taustd(df) = nanstd(goodph);
    [taumean(df), taustd(df), taubreakpt(df)] = circleMean(goodph, 1, []);
end

%add pts for filled area polygon for std deviations
yHi = taumean + taustd;
yLo = taumean - taustd;
yHi2 = taumean + 2*taustd;
yLo2 = taumean - 2*taustd;

fErrorBar = figure();
% pAvgSilico = errorbar(24*tau,taumean,taustd,'ks-',...
%     'linewidth',1,'markerfacecolor','k');
%hold on;
blueCol = cold(100);
lightBlue = [217 229 240]/255;
darkBlue = [191 209 229]/255;
pFill2 = fill(24*[tau fliplr(tau)],[yHi2 fliplr(yLo2)], lightBlue,...
    'linestyle','none');
hold on;
pFill = fill(24*[tau fliplr(tau)],[yHi fliplr(yLo)], darkBlue,...
    'linestyle','none','facealpha',0.5);
hold on;
pVitro=errorbar(24*inVitroDuty(4:end), inVitroDawnPhase(4:end), ...
    inVitroErr(4:end), 'bs',...
    'markersize',4,'markerfacecolor','b',...
    'markeredgecolor','none','linestyle','none');
hold on;
legend([pVitro,pFill, pFill2],...
    'in vitro','simulation \pm \sigma','simulation \pm 2\sigma',...
    'location','southwest');
legend boxoff;
xlabel('day length \tau (hours)');
ylabel('oscillator phase at dawn (rad/2\pi)');
set(gca,'xlim',[0 24],'ylim',[-0.5 0.5],'xtick',0:6:24,'ytick',[-0.6:0.2:0.6]);
grid off;
set(fErrorBar,'units','inches','position',[0 0 4 3]);
%%
if toexp_VitroSilico == 1
    export_fig([pwd '/' getDate('yyyy-mm-dd') '_fInVitroSilico_'...
        'NonParamBoot' num2str(max(boot)) '_ErrorBar_' ...
        getDate()],...
        '-cmyk','-painters','-pdf',fErrorBar);
end

%% proceed only if have a second simulation set to load (e.g., linearized step funs)
if TOLOAD_OLD_SIMULATION_TWO == 1
    SIMULATION_BOOT = load('saved_data/2017-03-30_09.56.39_NonParamBoot_Lin_I22_HiBrPt_1.mat');
    SIM = SIMULATION_BOOT.VAR;
    dawnphase = SIM.dawnphase;
    duskphase = SIM.duskphase;
    duskphaseshift = SIM.duskphaseshift;
    dawnphaseshift = SIM.dawnphaseshift;
    peakT = SIM.peakT;
    
    %% collect variables for plotting
    %plotphase(df,1:numel(numcyc),s)
    plotphase = dawnphase+dawnphaseshift;
    plotphase = mod(plotphase,1);
    plotphase(plotphase > 0.3) = plotphase(plotphase > 0.3) - 1;
    
    %% select only those tau's where entrained stably
    STDEV_LIM = 0.01;
    
    NUMENTR=4; %for nonlinear case
    for b=boot
        for df=1:numel(dutyfrac)
            plotphase_stdev(df,b) = std(squeeze(plotphase(df,(end-NUMENTR):end,b)));
            [~,plotphase_stdev(df,b),~] = ...
                circleMean(squeeze(plotphase(df,(end-NUMENTR):end,b)),1,[]);
            if plotphase_stdev(df,b) > STDEV_LIM
                entrained(df,b) = 0;
            else
                entrained(df,b) = 1;
            end
        end
    end
    
    %display some entrainment stats
    disp(['Total no. linear sims.: ' num2str(numel(entrained))]);
    disp(['No. entrained linear sims.: ' num2str(sum(sum(entrained)))]);
    disp(['Frac entrained linear: ' num2str(...
        sum(sum(entrained))/...
        numel(entrained)...
        )]);
    
    %specifically for tau between 8 and 16
    gf = dutyfrac <= 16/24 & dutyfrac >= 8/24;
    disp(['Total no. linear sims. tau=6-18: ' num2str(numel(entrained(gf,:)))]);
    disp(['No. entrained linear sims. tau=6-18: ' num2str(sum(sum(entrained(gf,:))))]);
    disp(['Frac entrained linear tau=6-18: ' num2str(...
        sum(sum(entrained(gf,:)))/...
        numel(entrained(gf,:))...
        )]);
    
    
    %% add in data to the first figure
    figure(fInVitroSilico); 
    for b=boot
        if sum(entrained(:,b) > 0) %& goodm(b) == 1
            entr = logical(entrained(:,b));
            pSilicoLin(b)=plot(24*dutyfrac(entr), plotphase(entr,end,b),...
                'gs-','markersize',4,...
                'markerfacecolor','g',...
                'linewidth',2);
            hold on;
            uistack(pSilicoLin(b),'bottom');
        end
    end
    legend([pVitroErr,pSilico(1),pSilicoLin(1)],...
        {'in vitro','full L & D sim.','linear L & D sim'},...
        'location','south');
    legend boxoff;
    
    if toexp_VitroSilico == 1
        export_fig([pwd '/' getDate('yyyy-mm-dd') '_fInVitroSilico_'...
            'NonParamBoot' num2str(max(boot)) '_ErrorBar_' ...
            getDate()],...
            '-cmyk','-painters','-pdf',fInVitroSilico);
    end
    
    %% plot avg phase+/-std for entrained tau's for second dataset
    %go through dfs, collect [tau, taustd]
    for df=1:numel(dutyfrac)
        tau(df) = dutyfrac(df);
        entr = logical(entrained(df,:));
        goodph = plotphase(df,end,entr);
        goodph = goodph(goodph < 0.33);
        %     taumean(df) = nanmean(goodph);
        %     taustd(df) = nanstd(goodph);
        [taumean(df), taustd(df), taubreakpt(df)] = circleMean(goodph, 1, []);
    end
    
    %go back to old figure
    figure(fErrorBar);
    %pErr2 = errorbar(24*dutyfrac,taumean,taustd,'gs-','linewidth',2);
    greenCol = [96 192 96]/255;
    pErr3=fill(24*[tau fliplr(tau)],[taumean+taustd fliplr(taumean-taustd)], ...
        greenCol, 'linestyle','none','facealpha',0.5);
    
    %get rid of 2sigma on nonparam
    delete(pFill2);
    
    legend([pVitro,pFill,pErr3],...%, pFill2],...
        'in vitro \pm \sigma','simulation \pm \sigma','linearized simulation \pm \sigma',...
        'location','southwest');
    uistack(pVitro,'top');
    legend boxoff;
    %%
    if toexp_VitroSilico == 1
        export_fig([pwd '/' getDate('yyyy-mm-dd') '_fInVitroSilicoAndLin_'...
            'NonParamBoot' num2str(max(boot)) '_ErrorBar_' ...
            getDate()],...
            '-cmyk','-painters','-pdf',fErrorBar);
    end
       
end
