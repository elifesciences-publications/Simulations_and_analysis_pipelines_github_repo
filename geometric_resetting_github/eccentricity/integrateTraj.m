function [tout,polarTraj,xyTraj] = integrateTraj(orbit, IC, tspan)
%Integrates dOmega/dt, dr/dt from initial condition IC=(r0,theta0) for time t
%tout = [t0; ... tend]; polarTraj = [r theta; ...]; xyTraj = [x y; ...];
%polarTraj computed relative to orbit. xyTraj is in true cartesian
%coordinates relative to the origin (0,0).

%empirically found that ode23tb works more reliably, even if more slowly. 

%For a=10, ecc=10:
%ode23tb = 0.67 min
%ode15s =  breaks
%ode23s = 1.2min

%For a=2, ecc=2:
%ode23tb = 0.41 min
%ode15s = 0.36 min
%ode23s = 0.62 min
tout = [];
polarTraj = [];

try
    [tout, polarTraj] = ode23tb(@(t, pt) ...
        ddt_ellipse_full(t, pt, orbit), tspan, IC, orbit);
catch ME
    warning(['ode23tb failed during integration. Will try with ode23s.']);
    try
        [tout, polarTraj] = ode23s(@(t, pt) ...
            ddt_ellipse_full(t, pt, orbit), tspan, IC, orbit);
    catch
        warning(['ode23s failed during integration too. Will try with ode15s.']);
        [tout, polarTraj] = ode15s(@(t, pt) ...
            ddt_ellipse_full(t, pt, orbit), tspan, IC, orbit);
    end
end
    
tout = tspan;
xyTraj = cartesianCoordOn(orbit,polarTraj);


end

