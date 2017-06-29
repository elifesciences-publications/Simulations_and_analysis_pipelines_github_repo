%2017-05-22, EL: modified for simulations with wD/wL = 0.94 (as in data)
%2017-03-29, EL: generate Fig. 4-figSup3
%dependencies: export_fig.m from FileExchange for figure export

TOEXP = 0; %export figure? (1=yes)
OMR = 0.93457;
%% compute m = d*(1-L)/(l+d*(1-L)) as a heatmap
Lrange = [0.0:0.02:1.0];
Drange = [0.01:0.02:1.01];
[Lslope, Dslope] = meshgrid(Lrange, Drange);
%Mslope = Dslope.*(1-Lslope)./(Lslope+Dslope.*(1-Lslope));
mslopefun = (@(L,D,OMR) (1- (-(1-L).*OMR./(D+L-L.*D)) - (1./(D+L-L.*D))));
Mslope = mslopefun(Lslope,Dslope,OMR);
Mslope(isinf(Mslope)) = nan;
Mslope(abs(Mslope) > 2) = nan;

%toss values less than 0 (set a certain threshold)
negLim = -0.00001;
Mslope(Mslope < negLim) = nan;

minM = min(min(Mslope));
maxM = max(max(Mslope));

Mslope(isnan(Mslope)) = min(min(Mslope))-0.03;

Ltick = Lrange(1):0.2:1; %Lrange(end);
Dtick = 0:0.2:1; %Drange(end);
%Dtick = Drange(1):0.2:Drange(end);

%slopefun
%mslopefun = (@(D,L) D*(1-L)/(L+D*(1-L)));
%mslopefun = (@(L,D,OMR) (1- (-(1-L)*OMR/(D+L-L*D)) - (1/(D+L-L*D))));

%isoclines for m = 0.4, 0.6
%dfun = (@(m,L) m.*L./((m-1).*(L-1.0)));
dfun = (@(m,L,OMR) (-OMR./(1.0-m)) - (1.0./((L-1.0).*(1.0-m))) + L./(L-1.0));
lgrid = 0:0.02:0.98;
d_lo = dfun(0.4,lgrid,OMR);
d_hi = dfun(0.6,lgrid,OMR);
d_vivo = dfun(0.53,lgrid,OMR);    %Fig. 2C in vivo curve
d_vitro = dfun(0.3762,lgrid,OMR); %Fig. 2C. avg. of fluor. m slopes

%box for our best guesses for l and d
% recall from step fun analysis
% from bootstrapping m along with l and d, for two values of OMR = wD/wL
%
% fitting over wide region of L and D
% -----------------------------------
% average boostrap estimate for all datasets
% boot lslope = 0.34244 +/- 0.027396
% boot dslope = 0.37518 +/- 0.048563
% boot mslope = 0.34451 +/- 0.042843, for OMR=0.93457

%indicate limits of crosshair (see printout of experimental_l_and_d_funs.m
%above)
nw = 'wide_fit_';
l_lims = 0.34244  + 0.027396*[-1 1];
d_lims = 0.37518 + 0.048563 *[-1 1];

%% plot figure
fM_LandD = figure();
colormap(jet);
imagesc(Lrange,Drange,Mslope);
hold on;

plot(lgrid,d_vivo,'linestyle','--','color',[0.7 0.7 0.7]);
hold on;
plot(lgrid,d_vitro,'linestyle','--','color',[0.7 0.7 0.7]);

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
set(gca,'xtick',Ltick,'ytick',Dtick,'xlim',[0 1], 'ylim',[0 1]);
set(gca,'ydir','normal');
axis('square');
xlabel('l');
ylabel('d');
grid off;
%box off;

colordata = colormap;
for i=1:1
    colordata(i,:) = [0.8 0.8 0.8];
end
colormap(colordata);

cbar=colorbar('location','westoutside');
hold on;
%plot(cbar,[0 1],0.5*[0 1],'w');
%set(cbar,'ylim',[0.0 1]);
set(cbar,'ytick',floor(10*[min(min(Mslope)):0.1:max(max(Mslope+0.2))])/10);
title(cbar,'m(l,d)');
grid off;

set(fM_LandD,'units','inches','position',[0 0 6 4]);

if TOEXP==1
    export_fig([getDate('yyyy-mm-dd') '_mSlopeHeatMap_' nw ...
        getDate('HH.MM.SS') '.pdf'],...
        '-cmyk','-painters',fM_LandD);
end