function [cartesianCoord] = cartesianCoordOn(orbit, polarCoord)
%Given polar coordinates relative to an orbit centered at (orbit.X,0) with
%radius orbit.R, convert them to cartesian coordinates relative to (0,0).

%if pass structure, return structure
if isstruct(polarCoord)
    x = orbit.X + (polarCoord.r).*cos(polarCoord.theta);
    y = (polarCoord.r).*sin(polarCoord.theta);
    
    cartesianCoord.x = x;
    cartesianCoord.y = y;

%if pass vectors, return vectors
%assume [r theta] column vecs
elseif isnumeric(polarCoord) 
    x = orbit.X + polarCoord(:,1).*cos(polarCoord(:,2));
    y = polarCoord(:,1).*sin(polarCoord(:,2));
    cartesianCoord = [x y];
else
    error(['cartesianCoordOn input of incorrect type!']);
end

end

