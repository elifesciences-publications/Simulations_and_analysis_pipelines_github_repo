%2016-09-18, EL: now plot fit and overlay normalized data from entire
%timeseries (even pts before fitting)
%2016-09-04, EL: modify to pass a vector of indices to resample data before
%plotting
%
% Dependencies: suplabel.m and export_fig.m from FileExchange to plot and
% save the figure.
%
function testParSet_Resampled(PARSET_hiATP, PARSET_loATP, SAMPLE_IND)
%test a parameter set by plotting
%   PARSET = parameter set, as 31x1 vector or {[1x4],...,[1x4]} cell array
%   PARSET_hiATP must contain fit params for sU1-sU9, sD10, in order.
%   PARSET_loATP must contain fit params for sD1-sD9, sU10, in order.

%% load expt data
INDIR = ['saved_data'];
INFILE = [INDIR '/2016-08-25_all_quants.xlsx'];

TOSAVE_EXAMPLE_FITFIG=0; %save figure with data + fit?

% need to trim non-fitting data points first, by analogy with
% bestFitNonParBoot.m

% script to load quants
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

%sample with replacement
AllRxns_trim = AllRxns_trim(SAMPLE_IND);
ActTime_trim = ActTime_trim(SAMPLE_IND);
pNet_trim = pNet_trim(SAMPLE_IND);

%% collect and arrange data

%put all high-ATP rxns into a single cell array
hiATPset = {'sU1','sU2','sU3','sU4','sU5','sU6','sU7','sU8','sU9','sD10'};
hiATPx = []; 
hiATPy = [];
hiATPstep = [];
[hiATPx, hiATPy, hiATPx_tofit, hiATPy_tofit, hiATPstep, hiATPlabel,...
    hiATPmu, hiATPsigma] = ...
    gatherToCell(hiATPset, stepUpTime, afterstep, ...
    AllRxns_trim, ActTime_trim, pNet_trim, TONORM);

%untrimmed, but normalized; afterstep=0, TONORM=0 here
[hiATPx_notr, hiATPy_notr, ~, ~, ~, ~, ~, ~] = ...
    gatherToCell(hiATPset, stepUpTime, 0, AllRxns, ActTime, pNet, 0);
for r=1:numel(hiATPy_notr)
    hiATPy_notr{r} = 100*(hiATPy_notr{r} - hiATPmu{r})./hiATPsigma{r};
end

%put all lo-ATP rxns into a single cell array
loATPset = {'sD1','sD2','sD3','sD4','sD5','sD6','sD7','sD8','sD9','sU10'};
loATPx = [];
loATPy = [];
[loATPx, loATPy, loATPx_tofit, loATPy_tofit, loATPstep, loATPlabel,...
    loATPmu, loATPsigma] = ...
    gatherToCell(loATPset, stepDownTime, afterstep, ...
    AllRxns_trim, ActTime_trim, pNet_trim, TONORM);

%untrimmed, but normalized; afterstep=0, TONORM=0 here
[loATPx_notr, loATPy_notr,~, ~, ~, ~, ~, ~] = ...
    gatherToCell(loATPset, stepDownTime, 0, AllRxns, ActTime, pNet, 0);

for r=1:numel(loATPy_notr)
    loATPy_notr{r} = 100*(loATPy_notr{r} - loATPmu{r})./loATPsigma{r};
end

%% check if PARSETs got passed as a {} or []
if isnumeric(PARSET_hiATP)
    fitPar_hiATP = convVecParToCA(PARSET_hiATP);
    fitPar_hiATP = fitPar_hiATP{1}; %kluge, but should only be one entry
else
    fitPar_hiATP = PARSET_hiATP;
end

if isnumeric(PARSET_loATP)
    fitPar_loATP = convVecParToCA(PARSET_loATP);
    fitPar_loATP = fitPar_loATP{1}; %kluge, but should only be one entry
else
    fitPar_loATP = PARSET_loATP;
end

%% go through fits to get finely sampled fit values for plotting
for r=1:numel(fitPar_hiATP)
    hiATPxx{r} = (min(hiATPx_tofit{r}):0.1:max(hiATPx_tofit{r}));
    hiATPyy{r} = sinusoidSimple(fitPar_hiATP{r}', hiATPxx(r), 0);
end

for r=1:numel(fitPar_loATP)
    loATPxx{r} = min(loATPx_tofit{r}):0.1:max(loATPx_tofit{r});
    loATPyy{r} = sinusoidSimple(fitPar_loATP{r}', loATPxx(r), 0);
end

%% plot stepUp rxns to check goodness of fit
fHighATP = figure();
for r=1:numel(fitPar_hiATP)-1
    axhiATP(r) = subplot(3,3,r);
    set(axhiATP(r),'xlim',[0 75],'ylim',[-10 80],...
        'xtick',0:12:84,'ytick',0:20:80);
    
    %control for stepUp (into hi-ATP rxns) is a rxn which is always in
    %loATP (never gets stepped up) = loATP{10}, which gets fit separately.
    %similarly, for step down the control is hiATP{10}, which got fit
    %separately. 
    pC = plot(loATPx{10},loATPy{10},'ks',...
        'markersize',4,'markerfacecolor','k','markeredgecolor','none');
    hold on;
    pCfit = plot(loATPxx{10},loATPyy{10},'k-','linewidth',2);
    
    
%     pSU(r) = plot(hiATPx{r},hiATPy{r},'bs',...
%         'markersize',4,'markerfacecolor','b','linewidth',2);
    pSUfit(r) = plot(hiATPxx{r},hiATPyy{r},'b-','linewidth',2);
    
    pSUnotr(r) = plot(hiATPx_notr{r}, hiATPy_notr{r},'bs',...
        'markersize',4,'markerfacecolor','b','markeredgecolor','none');
    
    if TONORM == 1
        set(axhiATP(r),'xlim',[0 75],'ylim',[-300 300],...
            'xtick',0:12:84,'ytick',-200:100:200);
    else
        set(axhiATP(r),'xlim',[0 75],'ylim',[-10 80],...
            'xtick',0:12:84,'ytick',0:20:80);
    end
    
    %plot light-dark boxes, after setting axes
    plotLDBoxes(axhiATP(r),0, hiATPstep{r},TONORM);
 
    %legend([pC pSU(r)],{'control',hiATPlabel{r}});
    title(hiATPlabel{r});
    grid off;
    legend boxoff;
end
suplabel('time (hours)','x');
suplabel('KaiC phosphorylation level (norm.)','y');
set(fHighATP,'units','inches','position',[0 0 10.5 7]); 

%% plot stepDown to check goodness of fit
floATP = figure();
for r=1:numel(fitPar_loATP)-1
    axloATP(r) = subplot(3,3,r);
    set(axloATP(r),'xlim',[0 75],'ylim',[-10 80],...
        'xtick',0:12:84,'ytick',0:20:80);
    pC = plot(hiATPx{10},hiATPy{10},'ks',...
        'markersize',4,'markerfacecolor','k','markeredgecolor','none');
    hold on;
    pCfit = plot(hiATPxx{10},hiATPyy{10},'k-','linewidth',2);
%     pSD(r) = plot(loATPx{r},loATPy{r},'rs',...
%         'markersize',4,'markerfacecolor','r','linewidth',2);
    pSDfit(r) = plot(loATPxx{r},loATPyy{r},'r-','linewidth',2);
    
    pSDnotr(r) = plot(loATPx_notr{r}, loATPy_notr{r},'rs',...
         'markersize',4,'markerfacecolor','r','markeredgecolor','none');
     
    if TONORM==1
        set(axloATP(r),'xlim',[0 75],'ylim',[-400 570],...
            'xtick',0:12:84,'ytick',-400:200:600);
    else
        set(axloATP(r),'xlim',[0 75],'ylim',[-10 80],...
            'xtick',0:12:84,'ytick',0:20:80);
    end
    
    %plot light-dark boxes, after axes lims have been set
    plotLDBoxes(axloATP(r),loATPstep{r},80,TONORM);
    
    %legend([pC pSD(r)],{'control',loATPlabel{r}});
    title(loATPlabel{r});
    grid off;
    legend boxoff;
end
suplabel('time (hours)','x');
suplabel('KaiC phosphorylation level (norm.)','y');
set(floATP,'units','inches','position',[0 0 10.5 7]); 

if TOSAVE_EXAMPLE_FITFIG==1
    export_fig([getDate('yyyy-mm-dd') '_stepUp_exBoot_' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',fHighATP);
    export_fig([getDate('yyyy-mm-dd') '_stepDown_exBoot_' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',floATP);
end



end

