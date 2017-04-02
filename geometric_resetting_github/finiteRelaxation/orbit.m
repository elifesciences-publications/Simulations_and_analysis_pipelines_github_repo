function [ output_args ] = untitled2( input_args )
%outline of entrainment simulation. no actual code.
%   Detailed explanation goes here


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate entrainment outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [t, [r theta], [x y], daynight] = entrain(dayOrb, nightOrb, tau, T, N)
%   --> [t, [r theta], [x y]] = integrateTraj(orbit, IC, t)
%           integrates dOmega/dt, dr/dt from initial condition IC for time t
%       --> [drdt, dthetadt] = ddt(r, theta, orbit)
%               compute ddt at polarCoord given the orbit that you're on
%               dtheta/dt = orbit.omega;
%               dr/dt = -orbit.a*(r-orbit.R);
% [orbit] = makeOrbit(R, X, a)
%   orbit.R = R; orbit.X = X; orbit.a = a;
% [r, theta] = polarCoordOn(orbit, [x y])
%   r = sqrt((x-orbit.X)^2 + y^2));
%   theta = arctan(y/(x-X)); %deal with branchcut too
% [x, y] = cartesianCoordOn(orbit, [r theta])
%   x = orbit.X + r.*cos(theta);
%   y = r.*sin(theta);