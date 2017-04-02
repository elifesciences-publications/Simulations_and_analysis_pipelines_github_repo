function [t, polarTraj, xyTraj, daynight] = ...
    entrain(dayOrb, nightOrb, tau, T, N, IC)
%Simulate entrainment to a day-night cycle of day length tau, period T for
%N cycles. During the day, oorbit around dayOrb; at dusk, orbit changes to
%the nightOrb.
%Supply IC=(r0,theta0) on dayOrb.

t=[];
polarTraj = [];
xyTraj = [];
daynight = [];
dt=0.1;

for n=1:N
    %first cycle is special b/c have empty variables
    if n==1
        %day
        [tout, pol, xy] = integrateTraj(dayOrb, IC, 0:dt:tau);
        t(:,1) = tout';
        polarTraj = pol;
        xyTraj = xy;
        daynight = ones(size(tout'));
        tout = []; pol = []; xy = [];
        
        %night
        %recompute initial condition here!
        nightStart = polarCoordOn(nightOrb, xyTraj(end,:));
        [tout, pol, xy] = integrateTraj(nightOrb, nightStart,(tau+dt):dt:T);
        polarTraj = [polarTraj; pol];
        xyTraj = [xyTraj; xy];
        daynight = [daynight; zeros(size(tout'))];
        t = [t; tout'];
        tout = []; pol = []; xy = [];
    
    %days 2 to N
    else
        %day. recompute initial condition in polarCoord relative to dayOrbit
        dayStart = polarCoordOn(dayOrb, xyTraj(end,:));
        [tout, pol, xy] = ...
            integrateTraj(dayOrb, dayStart, (n-1)*T+(dt:dt:tau));
        polarTraj = [polarTraj; pol];
        xyTraj = [xyTraj; xy];
        t = [t; tout'];
        daynight = [daynight; ones(size(tout'))];
        tout = []; pol = []; xy = [];
        
        %night.  recompute initial condition in polarCoord relative to nightOrbit
        nightStart = polarCoordOn(nightOrb, xyTraj(end,:));
        [tout, pol, xy] = ...
            integrateTraj(nightOrb, nightStart,(n-1)*T+((tau+dt):dt:T));
        polarTraj = [polarTraj; pol];
        xyTraj = [xyTraj; xy];
        t = [t; tout'];
        daynight = [daynight; zeros(size(tout'))];
        tout = []; pol = []; xy = []; 
    end
end

daynight = logical(daynight);

end

