%2017-03-28, EL: make figures for Fig6-figSup3(B).
%compute the bias in allocation of resources between dusk vs. dawn genes
%for different values of m.

function areaLL

close all;
clear all;

TOSAVE = 0;
TODISP = 1;

%timestep to integrate area in LL for dusk vs dawn genes
dx=1; %was 0.02
xs=-48:dx:72;

%load fits of dawn and dusk genes from Vijayan et al. (2009)
vijFit=load('2016-08-23_VijayanDawnDuskRawFit.mat', 'dawnFit', 'duskFit');
dawnFit=vijFit.dawnFit;
duskFit=vijFit.duskFit;

%give duskfit amplitude of 0.5, offset of 0
duskFit{1}(1) = 0;
duskFit{1}(3) = 0.5;
ys=(sinusoidSimple(duskFit{:}',{xs},0)+0.5);

CENTER = 18; 
dawngeneshift=(24/(2*pi))*(dawnFit{1}(4)-duskFit{1}(4));
duskgeneshift=0; 

xdusk = xs+duskgeneshift;
xdawn = xs+dawngeneshift;
 
xduskLarge = xdusk > 0 & xdusk < 24;
xdawnLarge_lo = xdawn < 24 & xdawn > 0;

%% helper functions
    function [dpstart, dpend] = getDPtimes(m,tau)
%        dpcenter = CENTER+(tau-12)*(m-0.5);
%         dpstart = [dpcenter-(24-tau)*0.5];
%         dpend = [dpcenter+(24-tau)*0.5];
        dpcenter = CENTER-(tau-12)*(m-0.5);
        dpstart = [dpcenter-(24-tau)*0.5];
        dpend = [dpcenter+(24-tau)*0.5];
    end

    function [nightTimes] = getNightTimes(m,tau)
        [dpstart, dpend] = getDPtimes(m,tau);
        nightTimes=[-24+(dpstart:dpend); dpstart:dpend; 24+(dpstart:dpend)];
    end
    
    function plotBoxes(m, tau, baseval)
        [nightTimes] = getNightTimes(m,tau);
        %plot LD gray boxes
        for nr=1:numel(nightTimes(:,1))
            a(nr)=area(nightTimes(nr,:),1.1*ones(size(nightTimes(nr,:)))+baseval,...
                baseval,'facecolor',[0.75 0.75 0.75],'edgecolor','none');
            set(gca,'layer','top');
            ch=get(a(nr),'children');
            %set(ch,'FaceAlpha',0.5)
            hold on;
            uistack(a(nr),'top');
        end
    end

    function plotDawnDuskLines(m,tau,baseval)
        [nightTimes] = getNightTimes(m,tau);
        duskCol = [0.75 0.75 0.75];
        for nr=1:numel(nightTimes(:,1))
            %dusk='k:',dawn='k--'
            plot(ones(1,10)*nightTimes(nr,1),...
                ones(1,10).*linspace(baseval,1.1+baseval,10),...
                ':','color',duskCol,'linewidth',1);
            hold on;
            plot(ones(1,10)*nightTimes(nr,end),...
                ones(1,10).*linspace(baseval,1.1+baseval,10),...
                '--','color',duskCol,'linewidth',1);            
        end
    end

    function plotWaveforms(baseval)
        duskCol = [0.2 0.2 0.2];
        dawnCol = [255 140 0]./255;
        %plot waveforms
        p(1)=plot(xdusk(xduskLarge),ys(xduskLarge)+baseval,'-','linewidth',2, ...
            'color',duskCol,'DisplayName','dusk gene');
        hold on;
        p(2)=plot(xdawn(xdawnLarge_lo),ys(xdawnLarge_lo)+baseval,'-','linewidth',2,...
            'color',dawnCol, 'DisplayName', 'dawn gene');
%         p(3)=plot(xdawn(xdawnLarge_hi),ys(xdawnLarge_hi)+baseval,'linewidth',2,...
%             'color',[0.2 0.2 0.2], 'DisplayName', 'dawn gene');
        set(gca,'layer','top');
        
    end

    function plotHorLines(baseval)
        blackCol = 'k';
        plot(xdusk(xduskLarge),ones(size(xdusk(xduskLarge))).*baseval,'-',...
            'color',blackCol,'linewidth',1);
    end
    
    function [dawngene, duskgene] = areaInLL(m, tau)
        [dpstart, dpend] = getDPtimes(m,tau);
        
        if dpstart > 24
            warning(['DP start time > 24 hrs!']);
        end
                        
        duskgene.LL = 0; 
        duskgene.DD = 0;
        dawngene.LL = 0; 
        dawngene.DD = 0;
        dispif(TODISP,['m=' num2str(m) ', tau=' num2str(tau)]);
        dispif(TODISP,['dpstart=' num2str(dpstart) ', dpend=' num2str(dpend)]);
        
       
        for i=1:numel(xdusk(xduskLarge))
            thisx = xdusk(xduskLarge);
            thisy = ys(xduskLarge);
            
            if dpstart <= 0 && (thisx(i) > dpend && thisx(i) < dpstart+24)
                duskgene.LL = duskgene.LL + thisy(i)*dx;
                dispif(TODISP,['LL dusk x:' num2str(thisx(i)) ', dusk y: ' num2str(thisy(i),'%2.2f')]);
            elseif dpstart > 0 && dpend <= 24 && (thisx(i) < dpstart || thisx(i) > dpend)
                duskgene.LL = duskgene.LL + thisy(i)*dx;
                dispif(TODISP,['LL dusk x:' num2str(thisx(i)) ', dusk y: ' num2str(thisy(i),'%2.2f')]);
            elseif dpstart > 0 && (dpend > 24) && ...
                    ((thisx(i) > mod(dpend,24) && thisx(i) < dpstart) ...
                    || (thisx(i) > dpend))
                duskgene.LL = duskgene.LL + thisy(i)*dx;
                dispif(TODISP,['LL dusk x:' num2str(thisx(i)) ', dusk y: ' num2str(thisy(i),'%2.2f')]);
            else
                duskgene.DD = duskgene.DD + thisy(i)*dx;
                dispif(TODISP,['DD dusk x:' num2str(thisx(i)) ', dusk y: ' num2str(thisy(i),'%2.2f')]);
            end
         end
         
         for i=1:numel(xdawn(xdawnLarge_lo))
             thisx = xdawn(xdawnLarge_lo);
             thisy = ys(xdawnLarge_lo);
             if dpstart <= 0 && (thisx(i) > dpend && thisx(i) < dpstart+24) 
                 dawngene.LL = dawngene.LL + thisy(i)*dx;
                 dispif(TODISP,['LL dawn x:' num2str(thisx(i)) ', dawn y: ' num2str(thisy(i),'%2.2f')]);
             elseif dpstart > 0 && dpend <= 24 && (thisx(i) < dpstart || thisx(i) > dpend)
                 dawngene.LL = dawngene.LL + thisy(i)*dx;
                 dispif(TODISP,['LL dawn x:' num2str(thisx(i)) ', dawn y: ' num2str(thisy(i),'%2.2f')]);
             elseif dpstart > 0 && (dpend > 24) && ...
                     ((thisx(i) > mod(dpend,24) && thisx(i) < dpstart) ...
                     || (thisx(i) > dpend))
                dawngene.LL = dawngene.LL + thisy(i)*dx;
                dispif(TODISP,['LL dawn x:' num2str(thisx(i)) ', dawn y: ' num2str(thisy(i),'%2.2f')]);
             else
                dawngene.DD = dawngene.DD + thisy(i)*dx;
                dispif(TODISP,['DD dawn x:' num2str(thisx(i)) ', dawn y: ' num2str(thisy(i),'%2.2f')]);
             end
         end         
    end

%% main figure with dawn/dusk gene area in LL
fRatioLL = figure();
tau=1:24;
mval=[0 0.5 1.0]; %which m values to compare: m=0 (dawn-tracking), m=0.5 (midday-tracking), m=1 (dusk-tracking)
%make backgrounds of dawn vs dusk-favored regions
bg(1)=area(0:24,1.1*ones(size(0:24))-1.1,-1.1,'facecolor',[0.75 0.75 0.75],...
    'edgecolor','none');
hold on;
bg(2)=area(0:24,1.1*ones(size(0:24)),0,'facecolor',[255 140 0]./255,...
    'edgecolor','none');

for j=1:numel(mval)
    m=mval(j);
    for ti=1:numel(tau)
        [dawngene, duskgene] = areaInLL(m, tau(ti));
        fracdusk(ti,j) = (duskgene.LL-dawngene.LL)/(duskgene.LL+dawngene.LL);
        %fracdusk(ti,j) = (duskgene.LL-dawngene.LL);
        %fracdusk(ti,j) = (duskgene.LL-dawngene.LL);
        dispif(TODISP,tau(ti));
        dispif(TODISP,fracdusk(ti,j));
    end
end

colors = {'k','r','k'};
linesym = {':','-','--'};
for j=1:numel(fracdusk(1,:))
    plt(j)=plot(tau,fracdusk(:,j),linesym{j},'color',colors{j},'linewidth',3);
    set(gca,'Layer','top')
    hold on;
    lab{j}=['m=' num2str(mval(j))];
end
hold on;


set(gca,'xlim',[6 18],'ylim',[-1.1 1.1],'xtick',0:6:24,'ytick',-1:0.5:1,...
    'xticklabel',0:6:24,'yticklabel',[]);
grid off; 

set(fRatioLL,'units','inches','position',[0 0 2.25 2.25]); %was 1.5x2.25


end