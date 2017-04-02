%2016-11-07, EL: make Fig. 6F
%plot limit cycles from Phong et al. (2013, PNAS)

%relies on export_fig.m from FileExchange

clear all;
close all;
clc;

nightCol = [0.2 0.2 0.2];
dayCol = [255 215 0]./255;

%load data
phong25 = loadConnieData('25');
phong100 = loadConnieData('100');

%%
fLimCyc = figure();
p25=plot(phong25.S, phong25.T,'s-','color',nightCol,...
    'markersize',4,'markerfacecolor',nightCol,'markeredgecolor','none',...
    'linewidth',1);
hold on;
p100=plot(phong100.S, phong100.T, 's-','color',dayCol,...
     'markersize',4,'markerfacecolor',dayCol,'markeredgecolor','none',...
     'linewidth',1);
xlabel('%pS431');
ylabel('%pT431');
legend([p100, p25], {'100% ATP', '25% ATP'},'fontsize',12);
legend boxoff;

set(gca, 'xlim', [0 40], 'xtick',[0:10:40], ...
    'ylim', [0 40], 'ytick', [0:10:40],'color','none');
grid off;

set(fLimCyc,'units','inches','position',[0 0 3.5 3],'color','w');
export_fig([getDate('yyyy-mm-dd') '_PhongData_' getDate() '.pdf'],...
    '-cmyk','-painters',fLimCyc);