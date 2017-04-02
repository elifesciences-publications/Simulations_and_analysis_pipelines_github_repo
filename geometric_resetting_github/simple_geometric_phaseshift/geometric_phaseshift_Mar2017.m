%2017-03-28, EL: generates Fig. 6A-B. Fig. 6-fig. sup. 1
% dependencies: stepDL.m, stepLD.m
% uses arrow.m and export_fig.m from FileExchange

close all;
clear all;

%export
TOEXP = 0;

%what is the radius of the night circle? how far apart are the circles? 
R=10^[0.7];
X=10^[0.5];

%at what angles (in degrees) do you want to draw L and D steps?
%NOTE: do not set theta==0 here. 360 is OK. 
thetaSteps = [1 60:60:360]; %in degrees

%colors of day/night circles and arrows
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
deltaTheta = 89;  %degrees; MJR has 14
Lstep = zeros(length(1:deltaTheta:360),2);
Dstep = zeros(length(1:deltaTheta:360),2);

%% simulate D step function
fLD_only = figure();
% boxmin = min(-rL*1.1+xL,-rD*1.1+xD);
% boxmax = max(max(xL+rL*1.1,xD + rD*1.1), max(rL,rD)*1.1);

%set axis limits
minX = min(-rL*1.1+xL,-rD*1.1+xD);
minY = min(-rL*1.1, -rD*1.1);
boxmin = min(minX, minY);

maxX = max(rL*1.1+xL,rD*1.1+xD);
maxY = max(1.1*[rL, rD]);
boxmax = max(maxX, maxY);

axis([boxmin boxmax boxmin boxmax]); 

axis square
hold on
viscircles([xL yL], rL * 1, 'edgecolor', 'b');
%viscircles([xL yL], rL * .95, 'color','w');
viscircles([xD yD], rD * 1, 'edgecolor','r');
%viscircles([xD yD], rD * .95, 'color','w');
%viscircles([xL yL], rL * .95, 'color','w');
viscircles([xL yL], rL * 0.025,'edgecolor','k'); %draw center dots
viscircles([xD yD], rL * 0.025,'edgecolor','k');
hold off

%%
fLtoD_arrows = figure();
axis([boxmin boxmax boxmin boxmax]); 
axis equal
hold on
viscircles([xL yL], rL * 1, 'edgecolor', dayCol);
%viscircles([xL yL], rL * .95, 'color','w');
viscircles([xD yD], rD * 1, 'edgecolor',nightCol);
%viscircles([xD yD], rD * .95, 'color','w');

i = 1;
for theta=thetaSteps
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
       arrow([x0 y0], [xf yf], ...
           'Length', 10,'Width',0.5,...
           'baseangle',75,'tipangle',10,...
           'color',arrowCol,'linewidth',0.5);
    end
    
    %only increment if took a step
    i = i + 1;
end

%delete last row if didn't get to it b/c deltaTheta doesn't divide into 360
if Dstep(end,1) < Dstep(end-1,1)
    Dstep(end,:) = [];
end
hold off

grid off;
axis off;
set(gca,'color','w');
figSize=[0 0 2*3.25 2*2.1];
figName=['_R' num2str(log10(R)) '_X' num2str(log10(X))];
set(fLtoD_arrows,'units','inches','color','w','position',figSize);

if TOEXP == 1
export_fig([getDate('yyyy-mm-dd') figName '_arrowsLD_' getDate('HH.MM.SS') '.pdf'], '-cmyk',...
    '-painters',fLtoD_arrows);
end

%% repeat for the L step (arrows are off)
fDtoL_arrows = figure();
% %set axis limits
% minX = min(-rL*1.1+xL,-rD*1.1+xD);
% minY = min(-rL*1.1, -rD*1.1);
% boxmin = min(minX, minY);
% 
% maxX = max(rL*1.1+xL,rD*1.1+xD);
% maxY = max(1.1*[rL, rD]);
% boxmax = max(maxX, maxY);

axis([boxmin boxmax boxmin boxmax]); 
axis equal;

hold on
viscircles([xL yL], rL * 1, 'edgecolor', dayCol);
%viscircles([xL yL], rL * .95, 'edgecolor','w');
viscircles([xD yD], rD * 1, 'edgecolor', nightCol);
%viscircles([xD yD], rD * .95, 'edgecolor','w');
arrowCol = 'b';

i = 1;
for theta=thetaSteps
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
       arrow([x0 y0], [xf yf],...
           'Length', 10,'Width',0.5,...
           'baseangle',75,'tipangle',10,...
           'color',arrowCol,'linewidth',1);
       %plot(x0,y0,'co');
       i = i + 1;
    end
   
end

grid off;
axis off;
set(gca,'color','w');
set(fDtoL_arrows,'units','inches','color','w',...
    'position',figSize);

if TOEXP == 1
export_fig([getDate('yyyy-mm-dd') figName 'arrowsDL_' getDate('HH.MM.SS') '.pdf'], '-cmyk',...
    '-painters',fDtoL_arrows);
end

%delete last row if didn't get to it b/c deltaTheta doesn't divide into 360
if Lstep(end,1) < Lstep(end-1,1)
    Lstep(end,:) = [];
end
hold off

%% plot one example of how angles are computed for one D fun

arrowCol = 'r';

%set axis limits
minX = min(-rL*1.1+xL,-rD*1.1+xD);
minY = min(-rL*1.1, -rD*1.1);
boxmin = min(minX, minY);

maxX = max(rL*1.1+xL,rD*1.1+xD);
maxY = max(1.1*[rL, rD]);
boxmax = max(maxX, maxY);

axis([boxmin boxmax boxmin boxmax]); 
hold off

fExample_Shift = figure();
axis([boxmin boxmax boxmin boxmax]); 
axis square
hold on
viscircles([xL yL], rL * 1, 'edgecolor', dayCol,'linewidth',2);
viscircles([xD yD], rD * 1, 'edgecolor',nightCol,'linewidth',2);

for theta=75
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
       arrow([x0 y0], [xf yf], 'Length', 20, 'Width', 2, 'color',arrowCol);
    end
end

%plot angle on L circle
plot([xL xL+rL],[yL, yL],'k','linewidth',1);
plot([xL x0], [yL, y0], 'k','linewidth',1);

%plot angle on D circle
plot([xD xD+rD], [yD yD], 'k','linewidth',1);
plot([xD xf], [yD yf], 'k','linewidth',1);
set(gca,'color','w');

%plot center of L circle
%plot(xL,yL,'o','markersize',4,'markerfacecolor','b','markeredgecolor','none');

axis off;
grid off;

set(fExample_Shift,'units','inches','position',[0 0 2*3.25 2*2.1],'color','w');

if TOEXP == 1
export_fig([getDate('yyyy-mm-dd') '_exampleDeltaTheta_' ...
    getDate() '.pdf'],'-cmyk','-painters',fExample_Shift);
end