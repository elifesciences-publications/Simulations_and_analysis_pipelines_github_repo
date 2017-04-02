%SHOULD PROBABLY DO ALL OF THIS ANALYTICALLY%

function [peakTime, time2circle, ptOnOrbit, appPh0] = ...
    phaseOnCircle(orbit, xyStartPt)
%Given a cartesian xyStartPt (x,y), integrate oscillator dynamics until the 
%oscillator has reached a stable circular orbit. Return peakTime = first 
%time when oscillator phase on orbit is pi/2; ptOnOrbit = polar coord when
%oscillator first hits the orbit; appPh0 = apparent phase at end of 
%entrainment, appPh0=reachPhase-time2reachPhase*omega; time2circle = time
%from start point to time of hitting the orbit.

tend = 100;
done = 0;
IC = polarCoordOn(orbit,xyStartPt);
tnet = 0;

%iteratively increase integration time until you've reached circular orbit
while ~done
    %integrate the trajectory in time
    [tout, polarTraj, xyTraj] = integrateTraj(orbit, IC, 0:0.01:tend);
          
    %when did the traj. hit the circular orbit
    [time2orbit_temp, ind2orbit] = time2stableCircle(orbit, tout, xyTraj);
       
    %if reached orbit, compute how much time remains till theta=0.5
    if ~isempty(time2orbit_temp)
        assert(time2orbit_temp >= 0); % time > 0!
        
        ptOnOrbit = polarTraj(ind2orbit,:);
        disp(['ptOnOrbit: ' num2str(cartesianCoordOn(orbit,ptOnOrbit))]);
        reachPhase = mod(ptOnOrbit(2),2*pi);
        if orbit.omega > 0
            pkPhase = pi/2;
            if reachPhase < pkPhase
                extraTime = (pkPhase - reachPhase)/orbit.omega;
            else
                extraTime = ((2*pi + pkPhase) - reachPhase)/orbit.omega;
            end
        else
            %if omega < 0, traveling backwards, so pkPhase must be <
            %reachPhase
            pkPhase = pi/2;
            if reachPhase < pkPhase
                extraTime = ((pkPhase - 2*pi) - reachPhase)/orbit.omega;
            else
                extraTime = (pkPhase - reachPhase)/orbit.omega;
            end
        end
        
        assert(extraTime >= 0); %must be positive!
        
        done = 1;
               
        time2circle = tnet + time2orbit_temp;
        tnet = tnet + time2orbit_temp + extraTime;
        
        %also compute apparent phase at end of entrainment based
        %on reachPhase and the amount of time it took to get there
        %apparent_ph0 = ph_T - omega*T;
        
        %%%%%%%%%%%%%%%%%%%% NOTE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% THIS IS SILLY ON A CIRCLE. appPh0 == ph at end of entrainment!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        appPh0 = reachPhase - (tnet-extraTime)*orbit.omega;
        appPh0 = mod(appPh0, 2*pi);
        
    %if haven't reached orbit, integrate more from last pt
    else
        IC = polarTraj(end,:);
        tnet = tnet + max(tout);
    end
end

%times must be positive!
assert(tnet >= 0);
assert(time2circle >= 0);

peakTime = tnet;
disp(['peakTime=' num2str(peakTime)]);
disp(['startingFrom = ' num2str(xyStartPt)]);

%estimate phase on a given orbit
%   --> get stable orbit (as input)
%   --> integrate from input phase in time until are convinced that you're
%   on stable trajectory
%   --> find first time when have phase=0.5 on stable trajectory
%   --> report this time mod period

end
 
% %store trajectory
% if isempty(xyTraj2Orbit)
%     xyTraj2Orbit = xyTraj;
%     t2orbit = tout;
% else
%     xyTraj2Orbit = [xyTraj2Orbit; xyTraj];
%     t2orbit = [t2orbit; tout];
% end