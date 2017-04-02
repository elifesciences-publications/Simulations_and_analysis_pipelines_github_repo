function [toutSimple,polarTrajSimple,xyTrajSimple] = integrateTraj(orbit, IC, tspan)
%Integrates dOmega/dt, dr/dt from initial condition IC=(r0,theta0) for time t
%tout = [t0; ... tend]; polarTraj = [r theta; ...]; xyTraj = [x y; ...];
%polarTraj computed relative to orbit. 

% [tout, polarTraj] = ode45(@(t, pt) ddt(t, pt, orbit), tspan, IC, orbit);
% tout = tspan;
% xyTraj = cartesianCoordOn(orbit,polarTraj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for this simple case, replace by explicit solution:
%theta(t) = omega*(t-t0) (sign of omega + or - = CCW or CW)
%r(t) = R + (r0-R)*e^(-at)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp(['tspan=' num2str(tspan)]);

%make sense of rotation determined by sign of omega
theta = IC(2) + (orbit.omega.*(tspan - tspan(1)));

r = orbit.R + (IC(1)-orbit.R).*exp(-orbit.a.*(tspan - tspan(1)));
polarTrajSimple = [r' theta'];
toutSimple = tspan;
xyTrajSimple = cartesianCoordOn(orbit,polarTrajSimple);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%in manuscript:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%see simulate_entrainment_11Jan2017.m
%stepLD.m, stepDL.m
 
%angle spanned by day vs night

%day = 360 * tau / T;
%night = 360 * (T - tau) / T;

%day
% theta = mod(theta + day,360);

%L -> D
% theta = stepLD(theta,R,X);

%night
% theta = mod(theta + night,360);

%D -> L
% theta = stepDL(theta,R,X);

%definitions of orbits
%xL = 0; rL = 1;
%xD = X; rD = R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

