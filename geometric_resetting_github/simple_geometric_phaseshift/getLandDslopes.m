%2017-03-28, EL: Used to compute slopes of L and D in Fig. 6-figSup1. 
%Compute slopes of L(theta) and D(theta) step functions as
%a function of center-to-center distances between the cycles (XIN) and the
%radius of the night cycle (RIN).
%This function does not always find the best region for a linear fit. Use
%this with manual supervision. 

%relies on fitToLine.m and arrow.m (from FileExchange).
function [LSLOPE, DSLOPE, LSTDEV, DSTDEV] = getLandDslopes(XIN, RIN) 
%use R=0.75, X=1.4 for data in Fig. 6

INDIR = '.';
cd(INDIR);

TOPLOT_LD = 1; %plot L and D functions with slopes of their linear regions highlighted?
TOPLOT = 1; %plot day and night limit cycles with arrows transitions?

R=RIN;
X=XIN;

nightCol = [0.2 0.2 0.2];
dayCol = [255 215 0]./255;
arrowCol = 'r';

% day limit cycle is centered at the origin and has unit radius
xL = 0;
yL = 0;
rL = 1;
% WLOG, night limit cycle has radius R and is displaced along the x-axis X
xD = -X;
yD = 0;
rD = R;

% simulation parameters
deltaTheta = 28;  %degrees; MJR has 14
Lstep = zeros(length(1:deltaTheta:360),2);
Dstep = zeros(length(1:deltaTheta:360),2);

%get axis limits
minX = min(-rL*1.1+xL,-rD*1.1+xD);
minY = min(-rL*1.1, -rD*1.1);
boxmin = min(minX, minY);

maxX = max(rL*1.1+xL,rD*1.1+xD);
maxY = max(1.1*[rL, rD]);
boxmax = max(maxX, maxY);


%% plot L to D transitions with arrows
if TOPLOT == 1
    fLtoD_arrows = figure();
    axis([boxmin boxmax boxmin boxmax]);
    axis square
    hold on
    viscircles([xL yL], rL * 1, 'edgecolor', dayCol);
    %viscircles([xL yL], rL * .95, 'color','w');
    viscircles([xD yD], rD * 1, 'edgecolor',nightCol);
    %viscircles([xD yD], rD * .95, 'color','w');
end

i = 1;
for theta=0:deltaTheta:359+deltaTheta
    % calculate initial point on day limit cycle
    x0 = xL + rL * cos(theta*2*pi/360);
    y0 = yL + rL * sin(theta*2*pi/360);
    
    % assuming instantaneous attraction to dark limit cycle, find new angle
    xjump = x0 - xD; %center of D circle is (xD, 0)
    yjump = y0;
    
    if abs(yjump) > 0
       newTheta = atand(yjump/xjump);
       xf1 = xD + rD * cos(newTheta*2*pi/360);
       yf1 = yD + rD * sin(newTheta*2*pi/360);
       xf2 = xD - rD * cos(newTheta*2*pi/360);
       yf2 = yD - rD * sin(newTheta*2*pi/360); 
       % i'm too lazy to figure out a good way to decide which branch of
       % arctangent should be used, so just test which is closest to the
       % original point.
       dist1 = (xf1 - x0)^2 + (yf1 - y0)^2;
       dist2 = (xf2 - x0)^2 + (yf2 - y0)^2;
       
       if dist1 < dist2
          xf = xf1;
          yf = yf1;
       else
          xf = xf2;
          yf = yf2;
          if newTheta > 0
             newTheta = newTheta + 180;
          else
             newTheta = newTheta - 180;
          end
          
       end
       
       Dstep(i,1) = theta;
       Dstep(i,2) = newTheta - theta;
       if Dstep(i,2) < -180
           Dstep(i,2) = Dstep(i,2) + 360;
       end
        if Dstep(i,2) > 180
           Dstep(i,2) = Dstep(i,2) - 360;
       end
      
       %plot jump
       %line([x0 xf], [y0 yf]);
       if TOPLOT == 1
           arrow([x0 y0], [xf yf], 'Length', 7,'Width',1,'color',arrowCol);
       end
    end
    
    %only increment if took a step
    i = i + 1;
    
end

%delete last row if didn't get to it b/c deltaTheta doesn't divide into 360
if Dstep(end,1) < Dstep(end-1,1)
    Dstep(end,:) = [];
end

if TOPLOT == 1
hold off
grid off;
axis off;
set(gca,'color','none');
set(fLtoD_arrows,'units','inches','color','w','position',[0 0 2*3.25 2*2.1]);
% export_fig([getDate('yyyy-mm-dd') '_arrowsLD_' getDate() '.pdf'], '-cmyk',...
%     '-painters',fLtoD_arrows);
end


%% plot D to L transitions with arrows

if TOPLOT == 1
arrowCol = 'b';

fDtoL_arrows = figure();
%set axis limits
minX = min(-rL*1.1+xL,-rD*1.1+xD);
minY = min(-rL*1.1, -rD*1.1);
boxmin = min(minX, minY);

maxX = max(rL*1.1+xL,rD*1.1+xD);
maxY = max(1.1*[rL, rD]);
boxmax = max(maxX, maxY);

axis([boxmin boxmax boxmin boxmax]); 
axis square;

hold on
viscircles([xL yL], rL * 1, 'edgecolor', dayCol);
%viscircles([xL yL], rL * .95, 'edgecolor','w');
viscircles([xD yD], rD * 1, 'edgecolor',nightCol);
%viscircles([xD yD], rD * .95, 'edgecolor','w');
end

i = 1;
for theta=0:deltaTheta:359
    % calculate initial point on night limit cycle
    x0 = xD + rD * cos(theta*2*pi/360);
    y0 = yD + rD * sin(theta*2*pi/360);
    
    % assuming instantaneous attraction to light limit cycle, find new angle
    xjump = x0; %center of L circle is (0,0)
    yjump = y0;
    
    if abs(yjump) > 0
       newTheta = atand(yjump/xjump);
       xf1 = xL + rL * cos(newTheta*2*pi/360);
       yf1 = yL + rL * sin(newTheta*2*pi/360);
       xf2 = xL - rL * cos(newTheta*2*pi/360);
       yf2 = yL - rL * sin(newTheta*2*pi/360); 
       % i'm too lazy to figure out a good way to decide which branch of
       % arctangent should be used, so just test which is closest to the
       % original point.
       dist1 = (xf1 - x0)^2 + (yf1 - y0)^2;
       dist2 = (xf2 - x0)^2 + (yf2 - y0)^2;
       
       if dist1 < dist2
          xf = xf1;
          yf = yf1;
       else
          xf = xf2;
          yf = yf2;
          if newTheta > 0
             newTheta = newTheta + 180;
          else
             newTheta = newTheta - 180;
          end
          
       end
       Lstep(i,1) = theta;
       Lstep(i,2) = newTheta - theta;
       if Lstep(i,2) <= -180
           Lstep(i,2) = Lstep(i,2) + 360;
       end
       if Lstep(i,2) > 180
           Lstep(i,2) = Lstep(i,2) - 360;
       end
       %plot jump
       if TOPLOT == 1
           arrow([x0 y0], [xf yf], 'Length', 7,'Width',1,'color',arrowCol);
       end
       i = i + 1;
    end
   
end

%delete last row if didn't get to it b/c deltaTheta doesn't divide into 360
if Lstep(end,1) < Lstep(end-1,1)
    Lstep(end,:) = [];
end

if TOPLOT == 1
hold off
grid off;
axis off;
set(gca,'color','none');
set(fDtoL_arrows,'units','inches','color','w','position',[0 0 2*3.25 2*2.1]);
% export_fig([getDate('yyyy-mm-dd') '_arrowsLD_' getDate() '.pdf'], '-cmyk',...
%     '-painters',fLtoD_arrows);
end


%% plot L and D in a 2pi periodic fashion

%find breakpt
[~,LbreakInd] = max(abs(diff(Lstep(:,2))));
[~,DbreakInd] = max(abs(diff(Dstep(:,2))));

[~,Lmid] = min(abs(diff(Lstep(:,2))));
[~,Dmid] = min(abs(diff(Dstep(:,2))));


MIDSEP = 2; %how far from this midpt to fit to line

Lstep(1:LbreakInd+1,3) = 1;
Lstep(LbreakInd+2:end,3) = 2;

Dstep(1:DbreakInd,3) = 1;
Dstep(DbreakInd+1:end,3) = 2;

Lstep2 = Lstep;
Lstep2(:,1) = Lstep2(:,1) + 360;
Lstep2(:,3) = Lstep2(:,3) + 1;
Lstep = [Lstep; Lstep2];

[~,ix] = sort(Lstep(:,1));
Lstep = Lstep(ix,:);
[Lreg, LregIx] = findPosSlopeReg(-Lstep(:,2));

Dstep2 = Dstep;
Dstep2(:,1) = Dstep2(:,1) + 360;
Dstep2(:,3) = Dstep2(:,3) + 1;
Dstep = [Dstep; Dstep2];

[~,ix] = sort(Dstep(:,1));
Dstep = Dstep(ix,:);
[Dreg, DregIx] = findPosSlopeReg(-Dstep(:,2));

%get L slope
step = Lstep;
step(:,1:2) = step(:,1:2) / 360;
step(:,2) = -step(:,2);

% FITLIM=0.15;
% fitreg = [step(:,2) < FITLIM & step(:,2) > -FITLIM & step(:,3) == 2];

% d2 = step(:,3) == 2;
% [~,Lmid] = min(abs(diff(Lstep(d2,2))));
% fitreg = find(step(:,3) == 2,1, 'first')+[Lmid-MIDSEP:1:Lmid+MIDSEP];
%fitreg = dL > 0 & step(:,3) == 1 || step(:,3) == 2;

fitreg = LregIx;
[coef,~]=fitToLine(step(fitreg,1),step(fitreg,2));
fitStep = polyval(coef,step(fitreg,1));
LSTDEV = std(step(fitreg,2)-fitStep);
LSLOPE = coef(1);
if TOPLOT_LD == 1
    fLfun = figure();
    for stepset=1:numel(unique(step(:,3)))
        thisset = step(:,3) == stepset;
        plotL(stepset)=plot(step(thisset,1), step(thisset,2),...
            'color','b','linewidth',1);
        hold on;
    end
    plot(step(fitreg,1), polyval(coef,step(fitreg,1)),'k--','linewidth', 2);
    title(['X=' num2str(XIN) '. R=' num2str(RIN) '. l=' num2str(coef(1),'%2.2f')]);
    xlabel('step phase (rad/2\pi)');
    ylabel('phase shift (rad/2\pi)');
    grid off;
    set(gca, 'xlim', [0 2], 'xtick', [0:0.5:2],...
        'ylim',[-1 0.5], 'ytick',[-0.5:0.25:0.5]);
    set(fLfun,'units','inches','position',[0 0 4 3],'color','w');
end

%get D slope
step = Dstep;
step(:,1:2) = step(:,1:2) / 360;
step(:,2) = -step(:,2);

% d2 = step(:,3) == 2;
% [~,Dmid] = min(abs(diff(Dstep(d2,2))));
% fitreg = [step(:,2) < FITLIM & step(:,2) > -FITLIM & step(:,3) == 2];
% fitreg = [Dmid-MIDSEP:1:Dmid+MIDSEP];

return;

fitreg = DregIx;
[coef,~]=fitToLine(step(fitreg,1),step(fitreg,2));
fitStep = polyval(coef,step(fitreg,1));
DSTDEV = std(step(fitreg,2)-fitStep);
DSLOPE = coef(1);
if TOPLOT_LD == 1
    fDfun = figure();
    for stepset=1:numel(unique(step(:,3)))
        thisset = step(:,3) == stepset;
        plotD(stepset)=plot(step(thisset,1), step(thisset,2),...
            'color','r','linewidth',1);
        hold on;
    end
    plot(step(fitreg,1), polyval(coef,step(fitreg,1)),'k--','linewidth', 2);
    title(['X=' num2str(XIN) '. R=' num2str(RIN) '. d=' num2str(coef(1),'%2.2f')]);
    grid off;
    set(gca, 'xlim', [0 2], 'xtick', [0:0.5:2],...
        'ylim',[-1 0.5], 'ytick',[-0.5:0.25:0.5]);
    xlabel('step phase (rad/2\pi)');
    ylabel('phase shift (rad/2\pi)');
    set(fDfun,'units','inches','position',[0 0 4 3],'color','w');
end

% export_fig([getDate('yyyy-mm-dd') '_LDfun_' ...
%     getDate() '.pdf'],'-cmyk','-painters',fLDfun);
