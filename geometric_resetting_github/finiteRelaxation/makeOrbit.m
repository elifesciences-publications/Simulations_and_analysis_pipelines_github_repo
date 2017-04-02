function [orbit] = makeOrbit(R, X, a, omega)
%Make orbit structure containing orbit radius R, center X, radial
%attraction term a.
%Assume orbit is centered on X axis, s.t. y coordinate is 0
   orbit.R = R; 
   orbit.X = X; 
   orbit.a = a;
   orbit.Y = 0;
   orbit.omega = omega;

end

