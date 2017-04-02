function [orbit] = makeEllipse(R, X, eccentricity, attraction, omega)
%Make structure for an elliptical containing orbit radius R, center X, 
%eccentricity delta and attraction towards orbit a. Freq. omega.
%Assume orbit is centered on X axis, s.t. y coordinate is 0
   orbit.R = R; 
   orbit.X = X; 
   orbit.a = attraction;
   orbit.Y = 0;
   orbit.omega = omega;
   orbit.delta = eccentricity;

   %eqn.:
   %x^2 + (delta^2)*(R^2)*y^2 = R^2
   
end

