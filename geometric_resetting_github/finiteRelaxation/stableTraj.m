%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REALIZED THIS WAS UNNECESSARY FOR CIRCULAR ORBITS %%%%%%%%
%%% USE stableCircleTraj(orbit) FUNCTION INSTEAD %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, polarTraj, xyTraj, time2stable] = stableTraj(orbit)
%Integrate orbit dynamics for a long time and return one full cycle on the
%orbit, in polar and cartesian coordinates. Also return approximate time it
%took to reach stable orbit from (0,0);

%try integrating for 10 cycles initially
time2stable = 240;
[tout, full_pol, full_xy] = integrateTraj(orbit, [0 0], 0:0.1:time2stable);
[yesno, meandist, xyStable, polarStable] = checkStable(2*pi/orbit.omega, tout,...
    full_xy, full_pol);

%consider adding a loop here to check how soon reach stability

%if don't reach stable orbit in 10 cycles, do 100
if yesno == 0
    time2stable = 2400;
    [tout, full_pol, full_xy] = integrateTraj(orbit, [0 0], 0:0.1:time2stable);
    [yesno, meandist, xyStable, polarStable] = checkStable(2*pi/orbit.omega, tout, ...
        full_xy, full_pol);
end

%package outputs
t = tout;
polarTraj = polarStable;
xyTraj = xyStable;

%plot to check
TOPLOT = 1;
if TOPLOT == 1
f=figure();
plot(full_xy(:,1), full_xy(:,2),'b');
hold on;
if ~isempty(xyStable)
    plot(xyStable(:,1), xyStable(:,2),'r');
end
title(['stable=' num2str(yesno) '. meandist = ' num2str(meandist,'%2.3f')]);
end

end

function [yesno, meandist, xyStable, polarStable] = ...
    checkStable(period, t, xyTraj, polarTraj)
%check if the xyTrajectory ([x1 y1; x2 y2; ...]) sampled at timepoints 
%given by t=[t1 t2 t3 ... ] is periodic with given period.
%Check this by computing Euclidean distance between the last two periods of
%xyTraj.
%Output boolean yesno (=1 if trajectories overlap, 0 otherwise) and mean
%distance between the two trajectories. Also return xy coords of stable
%traj in xyStable.

%need at least two full periods to compare
assert(max(t) - min(t) > 2*period);

%find last two periods
last = t > max(t) - period;
prelast = t <= max(t) - period & t > max(t) - 2*period;
if sum(last) ~= sum(prelast)
    if sum(last) > sum(prelast)
        numDiff = sum(last) - sum(prelast);
        ind = find(prelast == 1, 1, 'last');
        prelast(ind+1:ind+numDiff) = 1;
    else
        numDiff = sum(prelast) - sum(last);
        ind = find(last == 1, 1, 'last');
        last(ind+1:ind+numDiff) = 1;
    end
end

%both should be the same length
assert(sum(prelast) == sum(last));

%find xy trajectories from last and perlast periods
xyLast = xyTraj(last,:);
xyPrelast = xyTraj(prelast,:);

%compute Eucledian distances
xyDiff = xyLast - xyPrelast;
meandist = mean(sqrt(sum(xyDiff.^2,2)));

%output, threshold yesno
if meandist > 0.1
    yesno = 0;
    xyStable = [];
    polarStable = [];
else
    yesno = 1;
    xyStable = xyTraj(last,:);
    polarStable = polarTraj(last,:);
end

end