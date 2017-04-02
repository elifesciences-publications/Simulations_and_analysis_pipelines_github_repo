function [t, polarTraj, xyTraj] = circleTraj(orbit)
%Return polar and xy coordinates of one full cycle on a circle orbit.
%Integrate orbit dynamics explicitly.
%Return [t1; t2; ...], [r0 theta0; ...], [x0 y0; x1 y1; ...].

period = abs(2*pi/orbit.omega);

%start from r=R, theta = 0
xyStart = cartesianCoordOn(orbit, [0 0]);

%integrate
[t, polarTraj, xyTraj] = integrateTraj(orbit, [orbit.R 0], 0:0.1:period);

end

