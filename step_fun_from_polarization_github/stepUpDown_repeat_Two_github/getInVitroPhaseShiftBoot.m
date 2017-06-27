%2017-02-14, EL: modify for polarization dataset, don't count relative to
%min of oscillation, but relative to ph0 of sine.
%2016-08-31, EL: modify to accept a step function as an argument
%2016-08-28, EL: load boostrapped L and D funs
function [ PHASESHIFT ] = ...
    getInVitroPhaseShiftBoot(PERTURBPHASE, STEPFUN)
%NOTE: PERTURBPHASE = [0,1], STEPFUN = [0, 2*pi]!, PHASESHIFT = [0,1];

funname = mfilename('fullpath');
%disp(['Getting PhaseShift from ' funname]);

step = STEPFUN;

% convert input radians to [0,1]; add pi/2 = 0.25*(2pi) to count relative to min. pt of
% oscillation
%step.phase = 0.25+(step.phase/(2*pi));
step.phase = step.phase/(2*pi);
%up.phase(up.phase > 1.5) = up.phase(up.phase > 1.5)-2.0;
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

