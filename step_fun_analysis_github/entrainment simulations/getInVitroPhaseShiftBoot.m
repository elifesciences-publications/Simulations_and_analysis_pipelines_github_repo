%interpolate input step functions to get phase shifts at specified phases
function [ PHASESHIFT ] = ...
    getInVitroPhaseShiftBoot(PERTURBPHASE, STEPFUN)
%NOTE: PERTURBPHASE = [0,1], STEPFUN = [0, 2*pi]!, PHASESHIFT = [0,1];

funname = mfilename('fullpath');
%disp(['Getting PhaseShift from ' funname]);

step = STEPFUN;

% convert input radians to [0,1]; ph0= min. pt of oscillation
% don't need to adjust phase of step for fluor. pol.-derived functions,
% because the experimental_L_and_D_funs.m script already takes care of
% that.
step.phase = step.phase/(2*pi);
step.phaseShift = step.phaseShift/(2*pi);

% set interpolation type
INTERPTYPE = 'linear';

% make sure input phase is in [0,1]
PERTURBPHASE = mod(PERTURBPHASE,1);

PHASESHIFT = interp1(step.phase, step.phaseShift, PERTURBPHASE, INTERPTYPE);

if isnan(PHASESHIFT)
    funname = mfilename('fullpath');
    warning(['linear interp error at ph=' num2str(PERTURBPHASE)]);
    
    warning('trying spline');
    PHASESHIFT = interp1(step.phase, step.phaseShift, PERTURBPHASE, 'spline');
%     
%     if isnan(PHASESHIFT)
%         disp('trying nearest');
%         PHASESHIFT = interp1(step.phase, step.phaseShift, PERTURBPHASE, 'nearest');
%     end
    
    if isnan(PHASESHIFT)
        error([funname ': phaseShift=NaN!']);
    end
end

% report phaseshift
PHASESHIFT = -PHASESHIFT;

% test by plotting
TOPLOT=0;
if TOPLOT==1
    fTest = figure();
    pStep=plot(step.phase,step.phaseShift,'b.-');
    hold on;
    pIn=plot(PERTURBPHASE,-PHASESHIFT','ko');
    %legend([pU,pD,pIn],{'stepUp','stepDown','input'});
end

end

