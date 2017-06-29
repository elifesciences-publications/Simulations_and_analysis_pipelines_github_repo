%2017-06-26, EL: make Fig. 2C, use updated in vivo fits
%2017-05-24, EL: make Fig. 2C by adding on in vivo data 
%2017-03-30, EL: generate Fig. 2-figSup2.
% Perform a global fit of entrainment datasets collected using the
% fluorescence polarization probe and overlay the resulting dayTime vs. day
% length curves on top of SDS-PAGE data (as in Fig. 2C).
%
% dependencies: export_fig.m from FileExchange to export figure.
%

clear all; close all;

INDIR = '.';
cd([INDIR]);

%export?
TOEXP = 1; %save fig?
TOSAVE_FIT_PARAMS = 0; %save fits?

%% load fluorescence polarization data
JANFILE = '2017-02-01_inVitroPP_LL_PlateReader.xlsx';
MARFILE = '2017-03-07_inVitroPP_PlateReader.xlsx';
JANSHEETS = {'Sheet2'};
MARSHEETS = {'Sheet5'};
dsJan = joinSheets(JANFILE,JANSHEETS);
dsMar = joinSheets(MARFILE,MARSHEETS);

%% format well names to append labels and normalize datasets to 0mean, 1std

%deal with Jan
ds = dsJan;
dsName = 'Jan';
for i=1:numel(ds.Properties.VarNames) 
    if i>3 & numel(ds.Properties.VarNames{i}) <= 4
        name = ds.Properties.VarNames{i};
        new_name = [name '_' dsName];
        ds.Properties.VarNames{i} = new_name;
    end
end
%normalize to (0,1)
[nr,nc]=size(ds);
polmeanstd = zeros(2,nc);
polmeanstd = mat2dataset(polmeanstd); 
polmeanstd.Properties.VarNames = ds.Properties.VarNames;
for i=4:(nc-1)
    y = double(ds(:,i));
    polmeanstd(1:2,i) = dataset([nanmean(y) nanstd(y)]');
    y = (y-nanmean(y))./nanstd(y);
    ds(:,i) = mat2dataset(y);
end
dsJan = ds;

%deal with March dataset
ds = dsMar;
dsName = 'Mar';
%normalize to (0,1)
[nr,nc]=size(ds);
polmeanstd = zeros(2,nc);
polmeanstd = mat2dataset(polmeanstd); 
polmeanstd.Properties.VarNames = ds.Properties.VarNames;
for i=4:(nc-1)
    y = double(ds(:,i));
    polmeanstd(1:2,i) = dataset([nanmean(y) nanstd(y)]');
    y = (y-nanmean(y))./nanstd(y);
    ds(:,i) = mat2dataset(y);
end
for i=1:numel(ds.Properties.VarNames) 
    if i>3 & numel(ds.Properties.VarNames{i}) <= 4
        name = ds.Properties.VarNames{i};
        new_name = [name '_' dsName];
        ds.Properties.VarNames{i} = new_name;
    end
end
dsMar = ds;

%% fix times relative to last dawn
dsJan.time = dsJan.time + 1;
dsMar.time = dsMar.time + 0.5;

%% which wells to fit
janWells = {'B10_Jan','C10_Jan','D10_Jan','E10_Jan','F10_Jan','H10_Jan'};
janPP = [6 8 10 12 14 24];
janGoodPlt = janPP < 24; %good for plotting
janGoodFit = janGoodPlt; %good for fitting

marWells = {'L15_Mar','L16_Mar','L17_Mar','L18_Mar','L19_Mar','L20_Mar','L21_Mar'};
marPP = [6 8 10 12 14 16 24]-0.5;
marPP(end) = 24;
marGoodPlt = [marPP < 24];
marGoodFit = [marPP > 6 & marPP < 18];

combWells = {janWells{:}, marWells{:}};
combPP = [janPP marPP];

%indices of jan vs mar in array
janInd = 1:numel(janPP);
marInd = numel(janPP)+(1:numel(marPP));

marGoodPlt = [logical(zeros(size(janPP))) marGoodPlt];
marGoodFit = [logical(zeros(size(janPP))) marGoodFit];

%% extract X,Y data, set up fit funs
for w=1:numel(combWells)
    x=[];
    expset = combWells{w}(end-3:end);
    switch expset
        case '_Jan'
            x(:,1) = dsJan.time;
            x(:,2) = dsJan.(combWells{w});
            %toss first two hours
            goodx = x(:,1) > 2;
            x = x(goodx,:);
            perInds=numel(combWells)*3+1;
        case '_Mar'
            x(:,1) = dsMar.time;
            x(:,2) = dsMar.(combWells{w});
            %toss first two hours
            goodx = x(:,1) > 2;
            x = x(goodx,:);
            perInds=numel(combWells)*3+2;
    end
    
    %period index changes based on expt
    cntset=3*(w-1);
    fitdata.(combWells{w}) = x;
    fitfun.(combWells{w}) = ...
                        (@(b,X) b(cntset+1) + b(cntset+2)*...
                            sin((2*pi*X(:,1)/b(perInds)) + b(cntset+3))...
                            );
end

%% fit
%set up regression coefficients
b = zeros(1,numel(combWells)+2);

%indices of various fit params
npars = 3*numel(combWells)+2;
os = 1:3:npars;
Amps = 2:3:npars;
phs = 3:3:npars;
perInds = npars-[1 0]; %last two entries

% initial guesses
b0 = zeros(size(b));
b0(os) = 0;
b0(Amps) = 1;
b0(phs) = 0;
b0(perInds) = 24;

%lower bounds
lb = -Inf*ones(size(b0));
lb(phs) = -pi; %phi 
lb(perInds) = 20;
lb(Amps) = 0; %Ampl.

%upper bounds
ub = Inf*ones(size(b0));
ub(phs) = pi; %phi
ub(perInds) = 28;
ub(Amps) = Inf; %Ampl.

%settings
options = optimoptions('lsqnonlin','ScaleProblem','Jacobian');
options.MaxFunEvals = 10000;

% fit
[b, resnorm, residual, exitflag, output ,lambda, jacobian] = ...
    lsqnonlin(@globalMultiLinFitFun,b0,lb,ub,...
    options,...
    fitdata, fitfun);

%% test fit quality by plotting
funnames = fieldnames(fitfun);
setnames = fieldnames(fitdata);
assert(numel(setnames) == numel(funnames));
numsets = numel(funnames);

for n=1:numsets
    thisData = fitdata.(setnames{n});
    thisFun = fitfun.(funnames{n});
          
    subplot(3,5,n);
    
    plot(thisData(:,1), thisData(:,2),'bs','markersize',2,'linewidth',0.5);
    hold on;
    
    %polarization data is sampled every 0.25 hrs
    if thisData(2,1) - thisData(1,1) <= 1
        yfit = thisFun(b,thisData(:,1));
        plot(thisData(:,1),yfit,'r--');
        
    %interpolate finely for %PKaiC (every 4 hr sampled)
    elseif thisData(2,1) - thisData(1,1) > 1
        plot(thisData(:,1), thisData(:,2),'bs','markersize',4,'linewidth',2);
        xxfine = [0:0.1:72]';
        yfitfine = thisFun(b,xxfine);
        plot(xxfine,yfitfine,'r--');
    end
    
    title(strjoin(strsplit(funnames{n},'_')));
    
    set(gca,'ylim',[-3 3],'ytick',-4:1:4,'xlim',[0 72],'xtick',0:12:72);
    
    if n > 2
        xlabel('time (hours)');
    end
    
    ylabel('norm. FP');
end

set(gcf,'units','inches','position',[0 0 22 8.5]);

%export fit fig
exp_fig = @(fH,fName) export_fig([getDate() '_' fName '_' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf', fH);
fncs = {@(y) '', @(y) exp_fig(y{1},y{2})};
expif = @(x,y) feval(fncs{x+1},y);

expif(TOEXP,{[getDate('yyyy-mm-dd') '_fits_eachExpOwnPer_' ...
    getDate('HH.MM.SS') '.pdf'],gcf});

%% get fit quality
confint = nlparci(b,residual,'jacobian',jacobian);
sigma = abs(confint(:,1) - mean(confint,2));

%get variables for saving
out.BETA_FIT = b;
out.CHISQ = resnorm;
out.RESID = residual;
out.CI = confint; %checked that confint agrees with sigma

%% make var names for storage

for w=1:numel(combWells)
    expset = combWells{w}(end-3:end);
    switch expset
        case '_Jan'
            perInd = npars-1;
        case '_Mar'
            perInd = npars;
    end
    pkT(w) = (0.5*pi - mod(b(phs(w)),2*pi))*(b(perInd)/(2*pi));
    
    if pkT(w) < b(perInd) 
        pkT(w) = pkT(w) + b(perInd);
    end
    pkT_Err(w) = sigma(phs(w))*(b(perInd)/(2*pi));
    duskPhase(w) = b(phs(w)) + combPP(w)*2*pi/b(perInd);
end

fitPhase = b(phs);
fitPhase_Err = sigma(phs);
dawnPhase = fitPhase;

%do fits on pkT
[jan_pkT_lincoef, jan_pkT_errcoef] = ...
    fitToLineYErr(combPP(janGoodFit),pkT(janGoodFit),pkT_Err(janGoodFit))
[mar_pkT_lincoef, mar_pkT_errcoef] = ...
    fitToLineYErr(combPP(marGoodFit),pkT(marGoodFit),pkT_Err(marGoodFit))
%jan=0.39+/-0.06; mar=0.36+/-0.04

%do fits on dawnphase
[jan_ph_lincoef, jan_ph_errcoef] = ...
    fitToLineYErr(combPP(janGoodFit),...
        dawnPhase(janGoodFit),fitPhase_Err(janGoodFit)');
[mar_ph_lincoef, mar_ph_errcoef] = ...
    fitToLineYErr(combPP(marGoodFit),...
        dawnPhase(marGoodFit),fitPhase_Err(marGoodFit)');

%% save
vitro.pp = combPP;
vitro.fitPhase = fitPhase;
vitro.fitPhase_Err = fitPhase_Err;
vitro.dawnPhase = dawnPhase;
vitro.pkT = pkT;
vitro.pkT_Err = pkT_Err;
vitro.fitPer = [b(npars-1:npars)];
vitro.fitPer_Err = sigma((npars-1:npars));

vitro.marGoodPlt = marGoodPlt;
vitro.marGoodFit = marGoodFit;
vitro.janGoodPlt = janGoodPlt;
vitro.janGoodFit = janGoodFit;

vitro.jan_pkT_lincoef = jan_pkT_lincoef;
vitro.jan_pkT_errcoef = jan_pkT_errcoef;
vitro.mar_pkT_lincoef = mar_pkT_lincoef;
vitro.mar_pkT_errcoef = mar_pkT_errcoef;

vitro.jan_ph_lincoef = jan_ph_lincoef;
vitro.jan_ph_errcoef = jan_ph_errcoef;
vitro.mar_ph_lincoef = mar_ph_lincoef;
vitro.mar_ph_errcoef = mar_ph_errcoef;
    
if TOSAVE_FIT_PARAMS == 1
    save([getDate('yyyy-mm-dd') '_inVitroPP_fits_eachExpOwnPer_' ...
        getDate('HH.MM.SS') '.mat'], 'vitro');
end

%% load in vivo data
fittype='parabolic';
%disp(['vivo fittype = ' fittype]);
switch fittype
    case 'parabolic'
        vivoFile = ...
                ['Figures_and_Source_Data/'...
                 'in_vivo_entrainment_scaling_PARABOLA-26June2017-pk1.csv'];
        vivocoef = [0.5273    8.4990]; %from 2017-06-27 parabolic analysis (EL)
end

infile = vivoFile;
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%[^\n\r]'; %6 columns of doubles
fileID = fopen(infile,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
vivo = [dataArray{1:end-1}];
vivo = vivo(:,1:3); %in vivo is in first 3 columns

%% make figure
    
%compare with gel data
vitroGel = load('2017-03-30_inVitroPP_SDSPAGE.mat');
SDSPAGE = vitroGel.vitro;

fPkT_Compare = figure();
plot(SDSPAGE.pp,polyval(SDSPAGE.pkT_lincoef,SDSPAGE.pp),'k-');
hold on;
pPKaiC = errorbar(SDSPAGE.pp,SDSPAGE.pkT,SDSPAGE.pkT_Err,...
    's', 'color','k', 'markersize',4,'markerfacecolor','k',...
        'markeredgecolor','none','linewidth',1);

janCol = [0 255 255]/255;
plot(combPP(janGoodFit),polyval(jan_pkT_lincoef,combPP(janGoodFit)),...
    '-','color',janCol);
pJan=errorbar(combPP(janGoodPlt),pkT(janGoodPlt),pkT_Err(janGoodPlt),...
    's', 'color',janCol, 'markersize',4,'markerfacecolor',janCol,...
        'markeredgecolor','none','linewidth',1);

marCol = [0 200 255]/255;
plot(combPP(marGoodFit),polyval(mar_pkT_lincoef,combPP(marGoodFit)),...
    '-','color',marCol);
pMar=errorbar(combPP(marGoodPlt),pkT(marGoodPlt),pkT_Err(marGoodPlt),...
    'o', 'color', marCol, 'markersize',4,'markerfacecolor',marCol,...
        'markeredgecolor','none','linewidth',1);

%add in vivo
goodvivo = vivo(:,1) >=8 & vivo(:,1) <= 20;
plot(vivo(goodvivo,1),polyval(vivocoef,vivo(goodvivo,1)),'r-');
pVivo=errorbar(vivo(:,1),vivo(:,2),vivo(:,3),...
    's','color','r','markersize',4,'markerfacecolor','r',...
    'markeredgecolor','none','linewidth',1);

xlabel('day length \tau (hours)');
ylabel(['peak time t_{pk}' ' (hours after dawn)']);
axis([0 24 0 24]);
set(gca,'fontsize',12,'color','w','XTick',0:6:24,'YTick',0:6:24);
grid off;

% legend([pPKaiC, pJan,pMar],...
%     {'%P KaiC (trial 1)', 'KaiB pol. (trial 2)', 'KaiB pol. (trial 3)'},...
%     'location','southeast','fontsize',10);
% legend boxoff;

set(fPkT_Compare,'units','inches','Position',[0 0 3.3 3.3],'color','w');

expif(TOEXP,{[getDate('yyyy-mm-dd') '_compEntr_eachExpOwnPer_' fittype '_' ...
    getDate('HH.MM.SS') '.pdf'],gcf});

%% export combined dataset into a .csv file
ds_names = {'trial','day_length','peak_time','peak_time_error'};
gelTable = [1*ones(size((SDSPAGE.pp)')) (SDSPAGE.pp)' ...
    (SDSPAGE.pkT)' (SDSPAGE.pkT_Err)'];
janTable = [2*ones(size(combPP(janGoodPlt)'))  combPP(janGoodPlt)' ...
    pkT(janGoodPlt)' pkT_Err(janGoodPlt)'];
marTable = [3*ones(size(combPP(marGoodPlt)'))  combPP(marGoodPlt)' ...
    pkT(marGoodPlt)' pkT_Err(marGoodPlt)'];
vivoTable = [4*ones(size(vivo(:,1))) vivo];
dsAll = mat2dataset([gelTable; janTable; marTable; vivoTable]);
dsAll.Properties.VarNames = ds_names;
export(dsAll,'FILE','Fig2C-in_vitro_and_in_vivo_entrainment-June2017.csv','Delimiter',',');

%% export combined phase data into a .csv file
ds_names = {'trial','day_length','phase_at_dawn_rad','phase_error_rad'};
janTable = [2*ones(size(combPP(janGoodPlt)'))  combPP(janGoodPlt)' ...
   dawnPhase(janGoodPlt)' fitPhase_Err(janGoodPlt)];
marTable = [3*ones(size(combPP(marGoodPlt)'))  combPP(marGoodPlt)' ...
    dawnPhase(marGoodPlt)' fitPhase_Err(marGoodPlt)];
dsAll = mat2dataset([janTable; marTable]);
dsAll.Properties.VarNames = ds_names;
%export(dsAll,'FILE','in_vitro_entrainment_phase_units.csv','Delimiter',',');
