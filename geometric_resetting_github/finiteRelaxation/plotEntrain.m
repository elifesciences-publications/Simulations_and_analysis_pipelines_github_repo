function [orb, entrTraj, axlim] = plotEntrain(SimPick, SimVal, plotTimes, TOPLOT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nightCol = [0.2 0.2 0.2];
dayCol = [255 215 0]./255;

aval = SimVal.aval;
Rval = SimVal.Rval;
Xval = SimVal.Xval;
tauval = SimVal.tauval;
OMSIGN = SimVal.OMSIGN;
NUMCYC = SimVal.NUMCYC;
theta0 = SimVal.theta0;
PEAKT = SimVal.PEAKT;

pickA = SimPick.pickA;
pickR = SimPick.pickR;
pickX = SimPick.pickX;
pickTau = SimPick.pickTau;
pickTheta = SimPick.pickTheta;

%%

dO = makeOrbit(1, 0, aval(pickA),OMSIGN*2*pi/24);
nO = makeOrbit(Rval(pickR), Xval(pickX), aval(pickA), OMSIGN*2*pi/24);

[t, polTr, xyTr, dn] = ...
    entrain(dO, nO, tauval(pickTau), 24, NUMCYC, [dO.R theta0(pickTheta)]);
[tLL, polLL, release] = ...
    integrateTraj(dO, polarCoordOn(dO,xyTr(end,:)), ...
    0:1:squeeze(PEAKT(pickA,pickR,pickX,1,pickTau)));

%compute stable circular trajectories
[~,~,dCircle] = circleTraj(dO);
[~,~,nCircle] = circleTraj(nO);

%plot entrainment
%[255 180 100]/255
if isempty(plotTimes)
    lastDay = dn & t > (NUMCYC-1)*24;
    lastNight = ~dn & t > (NUMCYC-1)*24;
else
    assert(numel(plotTimes) == 2 & plotTimes(1) > 0 & ...
        plotTimes(2) > plotTimes(1));
    lastDay = dn & t > plotTimes(1) & t <= plotTimes(2);
    lastNight = ~dn & t > plotTimes(1) & t <= plotTimes(2);
end

orb.day = dCircle;
orb.night = nCircle;

entrTraj.xy = xyTr(lastDay | lastNight,1:2);
entrTraj.dn = dn(lastDay | lastNight); 
entrTraj.t = t(lastDay | lastNight);
entrTraj.xyLL = release;
entrTraj.tLL = tLL; 
entrTraj.polTr = polTr(lastDay | lastNight,1:2);
entrTraj.polLL = polLL;

%compute axis limits
minX = 1.25*min(-dO.R+dO.X,-nO.R+nO.X);
minY = 1.1*min(-dO.R, -nO.R);
boxmin = min(minX, minY);

maxX = 1.1*max(dO.R+dO.X,nO.R+nO.X);
maxY = 1.1*max([dO.R, nO.R]);
boxmax = max(maxX, maxY);

minX = min([min(xyTr(lastDay | lastNight,1)),...
    min(dCircle(:,1)) min(nCircle(:,1))]);
maxX = max([max(xyTr(lastDay | lastNight,1)),...
    max(dCircle(:,1)) max(nCircle(:,1))]);
minY = min([min(xyTr(lastDay | lastNight,2)),...
    min(dCircle(:,2)) min(nCircle(:,2))]);
maxY = max([max(xyTr(lastDay | lastNight,2)),...
    max(dCircle(:,2)) max(nCircle(:,2))]);
dX = maxX - minX;
dY = maxY - minY;

axlim = [minX-0.1*dX maxX+0.1*dX ...
    minY-0.1*dY maxY+0.1*dY];

%plot
if TOPLOT == 1
    
    %plot stable orbits
    plot(dCircle(:,1),dCircle(:,2), '-', 'color', dayCol, 'linewidth',6);   
    hold on;
    plot(nCircle(:,1),nCircle(:,2), '-', 'color', nightCol, 'linewidth',6);
    
    %plot entrainment
    plot(xyTr(lastDay,1),xyTr(lastDay,2),'color','r','marker','.',...
        'linestyle','none');
    hold on;
    plot(xyTr(lastNight,1),xyTr(lastNight,2),'color','r','marker','.',...
        'linestyle','none');
    
    %plot LL
%     plot(release(:,1), release(:,2), 'b--','linewidth',2);
    %plot(release(end,1), release(end,2), 'bx');
    
    set(gca, 'color',[0.9 0.9 0.9]);
    %set(gcf, 'color',[0.9 0.9 0.9]);
    
    axis(axlim);
    axis equal;
    grid off;
    
    %mark centers
%     plot(mean(dCircle(:,1)), mean(dCircle(:,2)),'x','color','r');
%     plot(mean(nCircle(:,1)), mean(nCircle(:,2)),'x','color','r');
%     disp([mean(nCircle(:,1)) mean(nCircle(:,2))]);
    
    title(['\tau=' num2str(tauval(pickTau))]);
    
end

end


% %set axis limits
% minX = 1.25*min(-dO.R+dO.X,-nO.R+nO.X);
% minY = 1.1*min(-dO.R, -nO.R);
% boxmin = min(minX, minY);
% 
% maxX = 1.1*max(dO.R+dO.X,nO.R+nO.X);
% maxY = 1.1*max([dO.R, nO.R]);
% boxmax = max(maxX, maxY);
% 
% axis([minX maxX minY maxY]);
% axis equal;
% grid off;
