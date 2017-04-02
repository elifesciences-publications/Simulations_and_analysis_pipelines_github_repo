function [t,ellipseXY,d2Ell] = ellipseFminbnd(pt,orbit,varargin)
%2017-03-10, EL
%Use Newton's method to find the nearest location on ellipse to the point
%pt.
%Based on: http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf

%to display search path of root-finding algorithm?
TOTEST = 0;

%to print intermediate outputs of this function?
TODISP = 0;
TODISPCORNER = 0;

%parse solver choice input: bisect or newton
solverChoice = 'fminbnd'; %looks like Matlab's default works best!
if numel(varargin) ~= 0
    solverChoice = varargin{1}; %
    %dispif(TODISP,solverChoice);
    if ~strcmp(solverChoice,'bisect') && ...
            ~strcmp(solverChoice,'newton') && ...
            ~strcmp(solverChoice,'fminbnd') 
        solverChoice = 'newton';
        warning('Didn''t understand solverChoice input! Using bisect.');
    end
end

% initialize outputs; useful in case fun breaks
t = [];
ellipseXY = [];
d2Ell = [];

%parse input parameters
dispif(TODISP, ['r, theta: (' num2str(pt(1)) ', ' num2str(pt(2)) ')']);

assert(pt(1) >= 0); % r must be positive!
r = pt(1);
theta = pt(2);

%(x,y) coordinates of query point (in the frame of reference of ellipse!)
x = r*cos(theta);
y = r*sin(theta);

dispif(TODISP, ['xy pt: (' num2str(x) ', ' num2str(y) ')']); %relative to center of orbit
cartXY = cartesianCoordOn(orbit,[pt(1) pt(2)]); %relative to true origin
dispif(TODISP, ['abs xy pt: (' num2str(cartXY(1)) ', ' num2str(cartXY(2)) ')']);

%get signs to determine which quadrant to search in
sx = sign(x);
sy = sign(y);

%convert pts to first quadrant pts
x = abs(x);
y = abs(y);

%ellipse parameters
R = orbit.R;
delta = orbit.delta;
%disp(['R=' num2str(R) ', delta=' num2str(delta)]); 
assert(R <= (1/delta)); %check that the major axis is along y axis

%circle case
if R == 1/delta
    dispif(1,['circle case, not origin']);
    if x > 0 || y > 0 %not origin
        ellipseXY = R*[x y]./sqrt(x.^2 + y.^2);
        ellipseXY = ellipseXY .* [sx sy];
        d2Ell = sqrt((ellipseXY(1)-x*sx).^2 + (ellipseXY(2)-y*sy).^2);
        if abs(ellipseXY(1)) > abs(sx*x) || ...
                abs(ellipseXY(2)) > abs(sy*y)  %query is inside ellipse
            d2Ell = -d2Ell;
            disp('corrected sign of distance');
        end
        t = atan2(ellipseXY(2),ellipseXY(1));
        return;
    else %could be at the origin
        dispif(1,['circle case, origin']);
        ellX = R;
        ellY = 0;
        if sx == 0 %at the origin, sign(r) = 0!
            ellipseXY = [ellX ellY];
        else
            ellipseXY = [ellX ellY].*[sx sy];
        end
        d2Ell = (sx*x - ellipseXY(1)); %already negative if inside ellipse
        dispif(1,['xaxis d2Ell=' num2str(d2Ell)]);
        plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, [], x,y, sx, sy);
        t = atan2(ellipseXY(2),ellipseXY(1));
        return;
    end
end

%parse corner cases
HAVECORNERS = 1; %consider corner cases?
if HAVECORNERS == 1
    if y > 0 
        if x > 0 %first quadrant
            %normal case, we're in first quadrant
            %Can consider the subcase for when you're on ellipse
            dispif(TODISPCORNER,'normal case');
            if ((x/R)^2 + (y*delta)^2 - 1) == 0 %on ellipse already!
                dispif(TODISPCORNER,'on ellipse');
                ellipseXY = [sx*x sy*y];
                d2Ell = 0;
                plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, [], x,y, sx, sy);
                t=atan2(ellipseXY(2),ellipseXY(1));
                return;
            end
            
        else % x == 0, y > 0, y axis
            dispif(TODISPCORNER,'on y axis');
            ellX = 0;
            ellY = 1/delta;
            ellipseXY = [ellX ellY].*[sx sy];
            d2Ell = sy.*y - sy.*ellY;
            plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, [], x,y, sx, sy);
            t = atan2(ellipseXY(2),ellipseXY(1));
            return;
        end
    else %y == 0, on x axis
        if (x < (R^2 - (1/delta)^2)/R)
            dispif(TODISPCORNER,'on x axis, case 1');
            ellX = (R*R*x)/(R*R - (1/delta)^2);
            ellY = (1/delta)*sqrt(1-(x/R)^2);
            d2Ell = sqrt((ellX-x)^2 + ellY^2);
                       
            if sx == 0 %at the origin, sign(r) = 0!
                ellipseXY = [ellX ellY];
            else
                ellipseXY = [ellX ellY].*[sx sy];
            end
            
            if abs(ellX) > abs(sx*x)  %query is inside ellipse
                d2Ell = -d2Ell;
                disp('corrected sign of distance');
            end
            
            plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, [], x,y, sx, sy);
            t = atan2(ellipseXY(2),ellipseXY(1));
            return;
        else
            dispif(TODISPCORNER,'on x axis, case 2');
            ellX = R;
            ellY = 0;
            if sx == 0 %at the origin, sign(r) = 0!
                ellipseXY = [ellX ellY];
            else
                ellipseXY = [ellX ellY].*[sx sy];
            end
            d2Ell = (sx*x - ellipseXY(1)); %already negative if inside ellipse
            dispif(TODISPCORNER,['xaxis d2Ell=' num2str(d2Ell)]);
            dispif(TODISPCORNER, ['abs xy pt: (' num2str(cartXY(1)) ', ' num2str(cartXY(2)) ')']);
            plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, [], x,y, sx, sy);
            t = atan2(ellipseXY(2),ellipseXY(1));
            return;
        end
    end
end
      
%decide whether to bisect and return or proceed with newton (lower)
switch solverChoice
    case 'bisect'
        [t,tcnt] = bisectRoot(orbit,x,y);
        [ellipseXY, d2Ell] = computeEllDist(t, R, delta, x, y, sx, sy);
        plotEllipse(TOTEST, orbit, abs(atan2(y,x)), ellipseXY, d2Ell, tcnt, x,y, sx, sy);
        return;
    case 'fminbnd'
        distformula = (@(t) sqrt((x-R*cos(t))^2 + (y-sin(t)/delta)^2));
        t_sln = fminbnd(distformula,0,pi/2); 
        [ellipseXY, d2Ell] = computeEllDist(t_sln, R, delta, x, y, sx, sy);
        %now need to adjust t based on the quadrant!!!
        t = atan2(ellipseXY(2),ellipseXY(1));
        dispif(TODISPCORNER,['t in fminbnd: ' num2str(t*(180/pi)) ' deg.']);
        plotEllipse(TOTEST, orbit, [], ellipseXY, d2Ell, t_sln, x, y, sx, sy);
        return;
end


%set up distance fun and its derivative for numerical solver
dfun = (@(t) sin(t).*cos(t)*((1/delta)^2 - R^2) + ...
             x.*R.*sin(t) - y.*(1/delta).*cos(t));
dprime = (@(t) (cos(t).^2 - sin(t).^2).*((1/delta)^2 - R^2) + ...
            x.*R.*cos(t) + y.*(1/delta).*sin(t));
        

%should limit search to a single quadrant of the ellipse!        
        
%set up initial guess
t0 = abs(atan2(R*y, (1/delta)*x)); %restrict to first quadrant, use atan

dispif(TODISP,['t0=' num2str(t0)]);

%use newton's method to get nearest pt
[t,tcnt] = newtonRoot(dfun, dprime, t0);
[ellipseXY, d2Ell] = computeEllDist(t, R, delta, x, y, sx, sy);
plotEllipse(TOTEST, orbit, t0, ellipseXY, d2Ell, tcnt, x,y, sx, sy);

%set up initial guess -- try another starting condition= 90-t0
t02 = (pi/2) - t0;

%use newton's method to get nearest pt
[t2,tcnt2] = newtonRoot(dfun, dprime, t02);
[ellipseXY2, d2Ell2] = computeEllDist(t2, R, delta, x, y, sx, sy);
plotEllipse(TOTEST, orbit, t02, ellipseXY2, d2Ell2, tcnt2, x,y, sx, sy);

%outputs -- choose smaller distance solution
if abs(d2Ell2) < abs(d2Ell)    
    t=t2;
    ellipseXY=ellipseXY2;
    d2Ell=d2Ell2;
end

%throw a warning if didn't get a solution
if isnan(d2Ell) || isempty(d2Ell)
    warning('d2Ell undefined!');
end

end

function plotEllipse(TOTEST, orbit, t0, ellipseXY, d2Ell, tcnt, x, y, sx, sy)
if TOTEST == 1
    figTest=figure();
    ts = [0:0.1:2*pi/orbit.omega];
    ell.x = orbit.R.*cos(ts);
    ell.y = sin(ts)./orbit.delta;
    plot(ell.x,ell.y,'m--');
    hold on;
    if ~isnan(tcnt)
        plot(sx*orbit.R*cos(tcnt), sy*(1/orbit.delta)*sin(tcnt),'cs-');
    end
    plot(ellipseXY(1), ellipseXY(2), 'go');
    plot(sx*x,sy*y,'bx');
    dispif(TOTEST, tcnt);
    axis('equal');
    if ~isempty(t0)
        title(['t0= ' num2str(t0*(180/pi)) ', dist = ' num2str(d2Ell, '%2.2f')]);
    else
        title(['no newton. , dist - ' num2str(d2Ell, '%2.2f')]);
    end
end
end

function [ellipseXY, d2Ell] = computeEllDist(t, R, delta, x, y, sx, sy)
%output negative distance if you're inside the ellipse

ellipseXY = [sx*R*cos(t) sy*(1/delta)*sin(t)];
d2Ell = sqrt((sx*x-ellipseXY(1))^2 + (sy*y-ellipseXY(2))^2);
%(x/R)^2 +(delta*Y)^2 = 1 on ellipse. if < 1, you're inside the ellipse
ellCond = (sx*x/R)^2 + (delta*sy*y)^2;
if ellCond < 1
    d2Ell = -d2Ell;
end
end

function [t_sln,t_all] = newtonRoot(dfun, dprime, tguess)

%display trial pts?
TODISP = 1;

%initialize search
t0 = tguess;

%set up Newton's method 
%(based on: https://m.njit.edu/Undergraduate/Matlab/M111MATLAB2S08/)
dt = 1;
f=dfun(t0);
Tol = 0.1; %was 0.001
count = 0;
tcnt = [];

while (dt > Tol || abs(f)>Tol)  %note that dx and f need to be defined for this statement to proceed
    count = count + 1;
    fprime = dprime(t0);
    dispif(TODISP, t0);
    f=dfun(t0);
    tnew = t0 - (f/fprime);   % compute the new value of x
    dt=abs(t0-tnew);          % compute how much x has changed since last step
    t0 = tnew;
    tcnt(count) = t0;
end

%output args
t_sln = t0;
t_all = tcnt;

end

function [t_sln,t_all] = bisectRoot(orbit,y0,y1)
%NOT used: wrapper for the distance-2-ellipse function; (y0,y1) = (x,y)
%coordinates of querty pt. orbit is passed to give the ellipse parameters.
e0 = orbit.R;
e1 = 1/orbit.delta;

z0 = y0/e0;
z1 = y1/e1;
g = (z0)^2 + (z1)^2 - 1;
r0 = (e0/e1)^2;

%outputs sbar, need to conver to x, y -> t
[s_sln,s_all] = getBisectRoot(r0,z0,z1,g);

y0_all = r0.*y0./(s_all + r0); %x
y1_all = y1./(s_all + 1); %y

t_all = tan(y1_all./y0_all);
t_sln = t_all(end);

end

function [s_sln, s_all] = getBisectRoot(r0, z0, z1, g)
%robust bisection algorithm from geometrictools.com, David Eberly

n0 = r0*z0;
s0 = z1 - 1;

if g < 0
    s1 = 0;
else 
    s1 = sqrt(n0^2 + z1^2) - 1;
end

s(1) = 0;
maxIter = 10;
Tol = 0.01;
for i=1:maxIter
    s(i) = (s0 + s1)/2;
    if (s(end) == s0) || (s(end) == s1) || abs(s1 - s0) < Tol
        break;
    end
    
    ratio0 = n0/(s(end)+r0);
    ratio1 = z1/(s(end)+1);
    g = (ratio0)^2 + (ratio1)^2 - 1;
    
    %bisect
    if g > 0
        s0 = s(end);
    elseif g < 0
        s1 = s(end);
    else
        break;
    end
end

s_sln = s(end);
s_all = s;

end

