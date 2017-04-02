function [dpt] = ddt_varspeed(t, pt, orbit)
%Compute d/dt at polar (r,theta) = (pt(1), pt(2)) computed relative to orbit.
%Used by ode45(). 
%Encodes equations of motions producing oscillations on a limit cycle, but
%with speed varying throughout the cycle.

%note:
%f = odefun(t,y), for a scalar t and a column vector y, 
% must return a column vector f corresponding to f(t,y).

r = pt(1);
theta = pt(2);

dtheta = orbit.omega + orbit.eps*sin(orbit.omega*t);

%note: this can be solved analytically easily!

dr = -(orbit.a).*(r-orbit.R);

dpt = [dr; dtheta];

end

