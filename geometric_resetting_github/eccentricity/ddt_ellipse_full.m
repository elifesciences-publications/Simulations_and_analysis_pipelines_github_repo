function [dpt] = ddt_ellipse_full(t, pt, orbit)
%2017-03-10, EL
%Compute d/dt at polar (r,theta) = (pt(1), pt(2)) computed relative to orbit.
%Used by ode45(). 
%Compute equations for an elliptical trajectory, any (r,theta) pt is
%attracted towards the nearest point on the ellipse.

%note:
%f = odefun(t,y), for a scalar t and a column vector y, 
% must return a column vector f corresponding to f(t,y).

r = pt(1);
theta = pt(2);

R = orbit.R;
delta = orbit.delta;

%[t, ellipseXY, d2Ell] = newtonEllipse(pt,orbit,'fminbnd');
[t, ellipseXY, d2Ell] = ellipseFminbnd(pt,orbit,'fminbnd');

% r_e = polarCoordOn(orbit,ellipseXY);
% ell_rdot = (-R*delta*y + y*orbit.omega)/cos(theta); %contribution to go around ellipse, not circle

%compute rdot based on r*rdot = x*xdot + y*ydot
%disp(['t=' num2str(t)]);
assert(isnumeric(t) & numel(t) == 1);

ell_rdot = -(1/2)*(( (cos(theta)/R)^2 + (delta*sin(theta))^2)^(-3/2))*  ...
    (2*orbit.omega)*(cos(theta)*sin(theta))*(-1/(R^2) + delta^2);

dtheta = orbit.omega;
dr = -(orbit.a).*(d2Ell) + ell_rdot; %pull tow. nearest pt on ellipse
%dr = -(orbit.a).*(r-r_e(1)); %pull tow. nearest pt on ellipse
%dr = cos(theta)*dtheta; %??? try this

TODISP = 0;
dispif(TODISP, ['ellipseXY: (' num2str(ellipseXY(1)) ...
                 ', ' num2str(ellipseXY(2)) ')']);
dispif(TODISP, ['d2Ell = ' num2str(d2Ell) ', ell_rdot = ' num2str(ell_rdot)]);
dispif(TODISP, ['dr, dtheta: (' num2str(dr) ', ' num2str(dtheta) ')']);
dispif(TODISP, char(10));

dpt = [dr; dtheta];

end

