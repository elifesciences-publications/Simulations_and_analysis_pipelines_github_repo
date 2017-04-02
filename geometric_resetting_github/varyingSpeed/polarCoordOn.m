function [polarPt] = polarCoordOn(orbit, cartesianPt)
%Convert cartesian point (x,y) to polar coordinates relative to a specific 
%Orbit centered at (orbit.X,0) with radius orbit.R. 
%If cartesian input is a vector cartesianPt = [x1 y1; x2 y2; ...],
%return polarPt = [r1 theta1; r2 theta2; ...].
    
%if pass structure, return structure
if isstruct(cartesianPt)
    x = cartesianPt.x;
    y = cartesianPt.y;
    
    r = sqrt((x-orbit.X).^2 + y.^2);
    theta = atan2(y,(x-orbit.X)); %deal with branchcut?

    polarPt.r = r;
    polarPt.theta = theta;
    
%if pass vector, return vector
%assume cartesianPt = [x1 y1; x2 y2; ...]
%return polarPt = [r1 theta1; r2 theta2; ...]
elseif isnumeric(cartesianPt)
    x = cartesianPt(:,1);
    y = cartesianPt(:,2);
    
    r = sqrt((x-orbit.X).^2 + y.^2);
    theta = atan2(y,(x-orbit.X)); %deal with branchcut?

    polarPt(:,1) = r;
    polarPt(:,2) = theta;
else
    error(['polarCoordOn input of incorrect type!']);
end


end

