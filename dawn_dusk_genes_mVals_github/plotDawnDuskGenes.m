%2017-03-28, EL: make figures for Fig6-figSup3(A)

function plotDawnDuskGenes

close all;
clear all;

%value of m to use to make Fig. 6-figSup3(A) panels
m=0.5;

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

CENTER = 18; %was 14
dawngeneshift=(24/(2*pi))*(dawnFit{1}(4)-duskFit{1}(4));
duskgeneshift=0; %was 8

xdusk = xs+duskgeneshift;
xdawn = xs+dawngeneshift;
 
xduskLarge = xdusk > 0 & xdusk < 48;
xdawnLarge_lo = xdawn < 48 & xdawn > 0;

%%
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
    
%% make figures showing which parts of transcriptional program are 
% "covered up" by night in different day-lengths (tau)

%figure
figWaves=figure();

%% top plot
tau=18;
BASEVAL=2.7;

plotWaveforms(BASEVAL);
plotBoxes(m, tau,BASEVAL-0.05);
plotDawnDuskLines(m,tau,BASEVAL-0.05);
plotHorLines(BASEVAL-0.1);

%% middle plot
tau=12;
BASEVAL=1.4;

plotWaveforms(BASEVAL);
plotBoxes(m,tau,BASEVAL-0.05);
plotDawnDuskLines(m,tau,BASEVAL-0.05);
plotHorLines(BASEVAL-0.1);

%% bottom plot
tau=6;
BASEVAL=0.1;

plotWaveforms(BASEVAL);
plotBoxes(m, tau,BASEVAL-0.05);
plotDawnDuskLines(m,tau,BASEVAL-0.05);

% pVertLine=plot(ones(1,10)*30,...
%                 ones(1,10).*[0:9],...
%                 ':','color','r','linewidth',1);
% uistack(pVertLine,'bottom');
            
%%
FONTSIZE=12;
xlabel('time (CT hours)','fontsize',FONTSIZE);
ylabel('expression level','fontsize',FONTSIZE);
%ylabel(['m_\tau=' sprintf('%2.1f',1-m)]);
set(gca,'xtick',0:6:48,'ytick',[]);
axis([0 48 0 3.9]);
%title(['m_\tau=' sprintf('%2.1f',1-m)]);
%legend(p,'location','northeastoutside');
%daspect([7.5 1 1]); %was 7.5x1x1

%set(gca,'linewidth',1,'color','none');
set(gca,'linewidth',1,'color','none');
set(gca,'fontsize',FONTSIZE);
box on;
grid off;

set(figWaves,'units','inches','position',[0 0 5 2.4],'color','w'); %was 1.5x2.25

end