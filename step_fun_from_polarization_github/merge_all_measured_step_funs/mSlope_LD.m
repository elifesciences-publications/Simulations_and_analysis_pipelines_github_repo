%2017-03-29, EL: generate Fig. 4-figSup3
%dependencies: export_fig.m from FileExchange for figure export

TOEXP = 0; %export figure? (1=yes)

%% compute m = d*(1-L)/(l+d*(1-L)) as a heatmap
Lrange = [0:0.02:1.0];
Drange = [0:0.02:1.0];
[Lslope, Dslope] = meshgrid(Lrange, Drange);
Mslope = Dslope.*(1-Lslope)./(Lslope+Dslope.*(1-Lslope));
Mslope(isinf(Mslope)) = nan;
Mslope(abs(Mslope) > 2) = nan;
Mslope(isnan(Mslope)) = min(min(Mslope))-0.03;

Ltick = Lrange(1):0.2:Lrange(end);
Dtick = Drange(1):0.2:Drange(end);

%slopefun
mslopefun = (@(D,L) D*(1-L)/(L+D*(1-L)));

%isoclines for m = 0.4, 0.6
dfun = (@(m,L) m.*L./((m-1).*(L-1.0)));
lgrid = 0:0.02:0.98;
d_lo = dfun(0.4,lgrid);
d_hi = dfun(0.6,lgrid);

%box for our best guesses for l and d
% recall from step fun analysis
% from bootstrapping m along with l and d, for two values of OMR = wD/wL
% lslope = 0.37832 +/- 0.053336
% dslope = 0.34435 +/- 0.033016
% mslope = 0.36366 +/- 0.056944, for OMR=1
l_lims = 0.37832 + 0.053336*[-1 1]; 
d_lims = 0.34435 + 0.033016*[-1 1];

%% plot figure
fM_LandD = figure();
imagesc(Lrange,Drange,Mslope);
hold on;

plot(lgrid,d_lo,'linestyle','--','color',[0.7 0.7 0.7]);
hold on;
plot(lgrid,d_hi,'linestyle','--','color',[0.7 0.7 0.]);

%% plot crosshair

%horiz.
chCol = 'w';
plot(l_lims, [1 1]*mean(d_lims),'-','color',chCol);
plot(mean(l_lims)+diff(l_lims)*0.15*[1 -1], ...
    [1 1]*(d_lims(1)),'-','color',chCol);
plot(mean(l_lims)+diff(l_lims)*0.15*[1 -1], ...
    [1 1]*(d_lims(2)),'-','color',chCol);

%vert
plot([1 1]*mean(l_lims), d_lims,'-','color',chCol);
plot([1 1]*(l_lims(1)), ...
    mean(d_lims)+diff(d_lims)*0.15*[1 -1], '-','color',chCol);
plot([1 1]*(l_lims(2)), ...
    mean(d_lims)+diff(d_lims)*0.15*[1 -1], '-','color',chCol);

 
ax(1)=gca;
set(gca,'xtick',Ltick,'ytick',Dtick);
set(gca,'ydir','normal');
axis('square');
xlabel('l');
ylabel('d');
grid off;
%box off;

colordata = colormap;
colordata(1,:) = [0.8 0.8 0.8];
colormap(colordata);

cbar=colorbar('location','westoutside');
hold on;
%plot(cbar,[0 1],0.5*[0 1],'w');
set(cbar,'ylim',[0.0 1]);
set(cbar,'ytick',floor(10*[min(min(Mslope)):0.1:max(max(Mslope+0.2))])/10);
title(cbar,'m(l,d)');
grid off;

set(fM_LandD,'units','inches','position',[0 0 6 4]);

if TOEXP==1
    export_fig([getDate('yyyy-mm-dd') '_mSlopeHeatMap_' ...
        getDate('HH.MM.SS') '.pdf'],...
        '-cmyk','-painters',fM_LandD);
end