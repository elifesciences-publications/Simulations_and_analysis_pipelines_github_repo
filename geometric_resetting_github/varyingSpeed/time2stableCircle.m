function [time2orbit, ind2orbit] = time2stableCircle(orbit, t, xyTraj)
%Given a circular orbit and timed xyTrajectory, return the time it takes for the
%trajectory to reach the stable orbit. 
%Return first timepoint when you've reached the orbit and the first index
%in timeseries/xyTraj when you're on the orbit. Empty if haven't reached
%the orbit.

%convert to polar
polarTraj = polarCoordOn(orbit, xyTraj);

%find when polar r coord. is close to orbit.R
ind2orbit = find(abs(polarTraj(:,1) - orbit.R) < 0.01,1,'first');

%return time too
if isempty(ind2orbit)
    time2orbit = [];
else
    time2orbit = t(ind2orbit);
end

end

