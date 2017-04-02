%2017-03-31, EL: generate Fig. 2-figSup1.
%2017-03-25, EL: modify to make a 2x2 supplemental figure. Fig2-fs1.
%2017-02-22, EL: global fit to polarization and %P KaiC data from
%out-of-phase reactions.
%
% dependencies: export_fig.m from FileExchange for figure export.
%

clear all; close all;

%export?
TOEXPORT = 0;

%% load %PKaiC data
%load all combined data %P KaiC SDS-PAGE data
nums = load(['2017-02-24_combinedGelNums_2017-02-24_17.56.45.txt']);
names = load(['2017-02-24_combinedGelNames2017-02-24_17.56.45.mat']);
names = names.final_names;
pkaic = 100-nums(:,4);
t1 = 20+[0 4 8 12 16 24];
t1 = [21.3075 25.2822 29.4172 33.3267 37.1425 46.1592]; %based on timestamps from the Excel files

gelt = [t1 t1 t1 t1 t1 t1]';

%normalize to (0,1)
pkaic_unnorm = pkaic;
kaimu = nanmean(pkaic);
kaisig = nanstd(pkaic);
pkaic = (pkaic - kaimu)./kaisig;

%% load fluorescence polarization data
INDIR = '.';
INFILE = '2016-02-19_PKaiC_Pol_PlateReader.xlsx';
cd([INDIR]);
INSHEETS = {'Sheet5','Sheet6','Sheet7','Sheet8','Sheet9','Sheet10','Sheet11'};
ds = joinSheets(INFILE,INSHEETS);
[nr,nc]=size(ds);

%%
polmeanstd = zeros(2,nc);
polmeanstd = mat2dataset(polmeanstd); 
polmeanstd.Properties.VarNames = ds.Properties.VarNames;

%normalize to (0,1)
for i=4:(nc-1)
    y = double(ds(:,i));
    polmeanstd(1:2,i) = dataset([nanmean(y) nanstd(y)]');
    y = (y-nanmean(y))./nanstd(y);
    ds(:,i) = mat2dataset(y);
end

%% set up regression coefficients
b = zeros(1,13);

%% fit functions
fitfun.polA = (@(b,X)  b(2) + b(3)* sin( (2*pi*X(:,1)/b(1)) + b(4)  ));
fitfun.polB = (@(b,X)  b(5) + b(6)* sin( (2*pi*X(:,1)/b(1)) + b(7)  ));
fitfun.gelA = (@(b,X)  b(8) + b(9)* sin( (2*pi*X(:,1)/b(1)) + b(10) ));
fitfun.gelB = (@(b,X)  b(11)+ b(12)*sin( (2*pi*X(:,1)/b(1)) + b(10) +(b(7)-b(4)) ));

%% put fitting data into correct format
% data = [t1 ydata1; t2 ydata2; ...];
ATPlevel = 100;
polAwell = 'G9';
polBwell = 'I9';
gelA = 1:6;
gelB = 7:12;

%polarization data
polTime = ds.time > 0;
data.polA = [double(ds.time(polTime)) double(ds(polTime,polAwell))];
data.polB = [double(ds.time(polTime)) double(ds(polTime,polBwell))];

% %PKaiC gel data
data.gelA = [gelt(gelA) pkaic(gelA)];
data.gelB = [gelt(gelB) pkaic(gelB)];

% un-normalized %PKaiC data
rawdata.gelA = [gelt(gelA) pkaic_unnorm(gelA)];
rawdata.gelB = [gelt(gelB) pkaic_unnorm(gelB)];

% un-normalized polarization data
mupolA = double(polmeanstd(1,polAwell));
sigpolA = double(polmeanstd(2,polAwell));
rawdata.polA = [double(ds.time(polTime)) ...
    double(ds(polTime,polAwell))*sigpolA+mupolA];

mupolB = double(polmeanstd(1,polBwell));
sigpolB = double(polmeanstd(2,polBwell));
rawdata.polB = [double(ds.time(polTime)) ...
    double(ds(polTime,polBwell))*sigpolB+mupolB];


%% fit

% initial guesses
b0 = [24	0	1	0	0	1	0	0	1	0	0	1];
lb = -Inf*ones(size(b0));
lb([4 7 10]) = -pi; %phi 
lb(1) = 20;
lb([3 6 9 12]) = 0; %Ampl.
ub = Inf*ones(size(b0));
ub([4 7 10]) = pi; %phi
ub(1) = 28;
ub([3 6 9 12]) = Inf; %Ampl.
options = optimoptions('lsqnonlin','ScaleProblem','Jacobian');
options.MaxFunEvals = 10000;

% fit
[b, resnorm, residual, exitflag, output ,lambda, jacobian] = ...
    lsqnonlin(@globalMultiLinFitFun,b0,lb,ub,...
    options,...
    data, fitfun);

%% make a single figure showing %PKaiC and polarizaton data for same rxns
% Combine %P KaiC data on a single plot (filled/open circles for the out of
% phase reactions, align based on fits) underneath polarization data.

fCombFig=figure();
subplot(2,2,1);
%plot polA. fit
yfit = fitfun.polA(b,data.polA(:,1));
plot(rawdata.polA(:,1),yfit*sigpolA+mupolA,'k-');
hold on;
%plot polA data
plot(rawdata.polA(:,1),data.polA(:,2)*sigpolA+mupolA,...
    'ks','markerfacecolor','k',...
    'markersize',2,'linewidth',0.5);
hold on;
ylabel('fluorescence polarization (mP)');
xlabel('time (hours)');
grid off;
set(gca,'xlim',[0 72],'xtick',0:12:72,'ylim',[120 180]);       

%gelA
subplot(2,2,3);
%plot gelA fit
xxfine = [0:0.1:72]';
yfitfine = fitfun.gelA(b,xxfine);
plot(xxfine,yfitfine*kaisig+kaimu,'k-');
hold on;
%plot gelA data
plot(data.gelA(:,1),data.gelA(:,2)*kaisig+kaimu,'ks','markersize',6,...
    'markerfacecolor','k','linewidth',0.5);
hold on;
set(gca,'ylim',[-5 100],'ytick',0:20:100,'xlim',[0 72],'xtick',0:12:72);
ylabel('%P KaiC');
xlabel('time (hours)');
grid off;

%polB
polBcol = 'b';
subplot(2,2,2);
%plot polB. fit
yfit = fitfun.polB(b,data.polB(:,1));
plot(rawdata.polB(:,1),yfit*sigpolB+mupolB,'-','color',polBcol);
hold on;
%plot polB data
plot(rawdata.polB(:,1),data.polB(:,2)*sigpolB+mupolB,...
    's','markeredgecolor',polBcol,'markerfacecolor',polBcol,...
    'markersize',2,'linewidth',0.5);
hold on;
ylabel('fluorescence polarization (mP)');
xlabel('time (hours)');
grid off;
set(gca,'xlim',[0 72],'xtick',0:12:72,'ylim',[120 180]);     

subplot(2,2,4);
%plot gelB data, shifted on X axis to match gelA data
%plot(data.gelB(:,1),data.gelB(:,2),'gs');
% shiftB = data.gelB(:,1);
% shiftB = shiftB+(b(1)/(2*pi)).*(b(7)-b(4));
%don't phase offset
yfitfine = fitfun.gelB(b,xxfine);
%plot gelB fit
plot(xxfine,yfitfine*kaisig+kaimu,'-','color',polBcol);
hold on;
%plot getB data
plot(data.gelB(:,1),data.gelB(:,2)*kaisig+kaimu,'s','color',polBcol,...
    'markersize',6,...
    'markerface',polBcol,'linewidth',0.5);
ylabel('%P KaiC');
xlabel('time (hours)');

set(gca,'ylim',[-5 100],'ytick',0:20:100,'xlim',[0 72],'xtick',0:12:72);
set(fCombFig,'units','inches','position',[0 0 7 5]);

grid off;

expif(TOEXPORT,{[getDate('yyyy-mm-dd') '_combinedPlt_' ...
    getDate('HH.MM.SS') '.pdf'],gcf});

%% arrange data to save into a .csv file
dsPol = [rawdata.polA(:,1) data.polA(:,2)*sigpolA+mupolA ...
    rawdata.polB(:,1) data.polB(:,2)*sigpolB+mupolB];
dsPol_names = {'rxnA_time_hrs','rxnA_polarization_mP',...
    'rxnB_time_hrs','rxnB_polarization_mP'};
dsPol = mat2dataset(dsPol);
dsPol.Properties.VarNames = dsPol_names;
export(dsPol,'FILE','polarization_reporter_raw_data.csv','Delimiter',',');

dsGel = [data.gelA(:,1) data.gelA(:,2)*kaisig+kaimu ...
    data.gelB(:,1) data.gelB(:,2)*kaisig+kaimu];
dsGel_names = {'rxnA_time_hrs','rxnA_gel_percentKaiCphospho',...
    'rxnB_time_hrs','rxnB_gel_percentKaiCphospho'};
dsGel = mat2dataset(dsGel);
dsGel.Properties.VarNames = dsGel_names;
export(dsGel,'FILE','percent_phospho_KaiC_raw_data.csv','Delimiter',',');


