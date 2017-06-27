%2017-02-13, EL: modify for plate reader
%2016-09-18, EL: now plot fit and overlay normalized data from entire
%timeseries (even pts before fitting)
%2016-09-04, EL: modify to pass a vector of indices to resample data before
%plotting
function testParSet_Resampled(PARSET_hiATP, PARSET_loATP, SAMPLE_IND)
%test a parameter set by plotting fits on top of raw data

%% load raw data, step times; set normalization/fit parameters
INDIR = ['saved_data'];
INFILE = [INDIR '/2017-02-24_stepDown_try2_PlateReader.xlsx'];
INSHEETS = {'Sheet3'};

TORESAMPLE = 0; %if 0, don't bootstrap -- use entire dataset

TOSAVE_EXAMPLE_FITFIG=1; %save figure with best fits?

% script to load plate reader data
[AllRxns, AllWells, ActTime, PolData] = loadPlateReader(INFILE,INSHEETS);

%'stepTimes' sheet has stepTimes relative to t=0 on plate reader
INFILE_STEP = [INDIR '/2017-02-23_stepDown_staged_schedule.xlsx'];
[stepTimes,~,~] = xlsread(INFILE_STEP,'stepTimes','B12:B27');
stepUpTime = [nan; nan]; 
stepDownTime = [stepTimes([1:7 9:15])];

afterstep = 2; %how much time to wait after step to start fitting
TONORM = 1; %normalize to mean0,std100 or to keep as raw data for fits/plots
TONORMTRIM = 0; %don't normalize while trimming
%%
allRxnTypes = {'sD1','sD2','sD3','sD4','sD5','sD6','sD7','sD8','sD9',...
    'sD10','sD11','sD12','sD13','sD14','sD15','sD16'};
         
[ActTime_trim, PolData_trim, AllRxns_trim, AllWells_trim] = ...
    gatherToMat(allRxnTypes, 2, AllRxns, AllWells,ActTime, PolData, TONORMTRIM);
AllRxns_trim = AllRxns_trim'; %must transpose
AllWells_trim = AllWells_trim';

%sample with replacement
AllRxns_trim = AllRxns_trim(SAMPLE_IND);
ActTime_trim = ActTime_trim(SAMPLE_IND);
PolData_trim = PolData_trim(SAMPLE_IND);

%% collect and arrange data

%put all high-ATP rxns into a single cell array
hiATPset = {'sD8','sD16'};
hiATPx = []; 
hiATPy = [];
hiATPstep = [];
[hiATPx, hiATPy, hiATPx_tofit, hiATPy_tofit, hiATPstep, hiATPlabel,...
    hiATPmu, hiATPsigma] = ...
    gatherToCell_Two(hiATPset, stepUpTime, afterstep, ...
    AllRxns_trim, ActTime_trim, PolData_trim, TONORM);

%untrimmed, but normalized; afterstep=0, TONORM=0 here
[hiATPx_notr, hiATPy_notr, ~, ~, ~, ~, ~, ~] = ...
    gatherToCell_Two(hiATPset, nan(size(stepUpTime)), 0, ... %pass NaNs as stepTimes to collect all data
    AllRxns, ActTime, PolData, 0);
for r=1:numel(hiATPy_notr)
    hiATPy_notr{r} = 100*(hiATPy_notr{r} - hiATPmu{r})./hiATPsigma{r};
end

%put all lo-ATP rxns into a single cell array
loATPset = {'sD1','sD2', 'sD3', 'sD4', 'sD5', 'sD6', 'sD7',...
            'sD9','sD10','sD11','sD12','sD13','sD14','sD15'};
loATPx = [];
loATPy = [];
[loATPx, loATPy, loATPx_tofit, loATPy_tofit,loATPstep, loATPlabel, ...
    loATPmu, loATPsigma] = ...
    gatherToCell_Two(loATPset, stepDownTime, afterstep, ...
    AllRxns_trim, ActTime_trim, PolData_trim, TONORM);

%untrimmed, but normalized; afterstep=0, TONORM=0 here
[loATPx_notr, loATPy_notr,~, ~, ~, ~, ~, ~] = ...
    gatherToCell_Two(loATPset, nan(size(stepDownTime)), 0,... %pass NaNs as stepTimes to collect all data
    AllRxns, ActTime, PolData, 0);
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
for r=1:numel(fitPar_hiATP)
       
    axhiATP(r) = subplot(1,2,r);
    set(axhiATP(r),'xlim',[0 75],'ylim',[-10 80],...
        'xtick',0:12:84,'ytick',0:20:80);
    
    %control for stepUp (into hi-ATP rxns) is a rxn which is always in
    %loATP (never gets stepped up) = loATP{10}, which gets fit separately.
    %similarly, for step down the control is hiATP{10}, which got fit
    %separately.    
    
%     pSU(r) = plot(hiATPx{r},hiATPy{r},'bs',...
%         'markersize',4,'markerfacecolor','b','linewidth',2);

    pSUfit(r) = plot(hiATPxx{r},hiATPyy{r},'b-','linewidth',2);
    hold on;   
    pSUnotr(r) = plot(hiATPx_notr{r}, hiATPy_notr{r},'bs',...
        'markersize',6); %,'markerfacecolor','b','markeredgecolor','none');
    
    if TONORM == 1
        set(axhiATP(r),'xlim',[0 60],'ylim',[-300 300],...
            'xtick',0:12:84,'ytick',-200:100:200);
    else
        set(axhiATP(r),'xlim',[0 60],'ylim',[-10 80],...
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
for r=1:numel(fitPar_loATP)
    cix = floor(r/8)+1; %control is sD8 for 1-7, sD16 for 8-15.
    
    axloATP(r) = subplot(4,4,r);
    set(axloATP(r),'xlim',[0 75],'ylim',[-10 80],...
        'xtick',0:12:84,'ytick',0:20:80);
    pC = plot(hiATPx{cix},hiATPy{cix},'k.',...
        'markersize',6);
    hold on;
    pCfit = plot(hiATPxx{cix},hiATPyy{cix},'k-','linewidth',2);
%     pSD(r) = plot(loATPx{r},loATPy{r},'rs',...
%         'markersize',4,'markerfacecolor','r','linewidth',2);
    pSDfit(r) = plot(loATPxx{r},loATPyy{r},'r-','linewidth',2);
    
    pSDnotr(r) = plot(loATPx_notr{r}, loATPy_notr{r},'r.',...
         'markersize',6);%,'markerfacecolor','r','markeredgecolor','none');
     
    if TONORM==1
        set(axloATP(r),'xlim',[0 60],'ylim',[-300 300],...
            'xtick',0:12:84,'ytick',-400:200:600);
    else
        set(axloATP(r),'xlim',[0 60],'ylim',[-10 80],...
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
suplabel('polarization (norm.)','y');
set(floATP,'units','inches','position',[0 0 10.5 7]); 

if TOSAVE_EXAMPLE_FITFIG==1
    export_fig([getDate('yyyy-mm-dd') '_stepUp_' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',fHighATP);
    export_fig([getDate('yyyy-mm-dd') '_stepDown_' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',floATP);
end



end

