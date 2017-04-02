function [time2orbit, ind2orbit] = time2stableEllipse(orbit, t, xyTraj)
%Given an elliptic orbit orbit and timed xyTrajectory, return the time it takes for the
%trajectory to reach the stable orbit. 
%Return first timepoint when you've reached the orbit and the first index
%in timeseries/xyTraj when you're on the orbit. Empty if haven't reached
%the orbit.

%convert to polar
polarTraj = polarCoordOn(orbit, xyTraj);

for i=1:numel(polarTraj(:,1))
    [s,eXY,d2E] = ellipseFminbnd(polarTraj(i,:),orbit,'fminbnd');
    ssave{i} = s;
    ellXY(i,:) = eXY;
    d2Ell(i,:) = d2E;
end

%find when polar r coord. is close to orbit.R
ind2orbit = find(abs(d2Ell) < 0.01,1,'first');

%return time too
if isempty(ind2orbit)
    time2orbit = [];
else
    time2orbit = t(ind2orbit);
end

end

