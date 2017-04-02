%% Fit TopCount trajectories to sinusoids with same period
%varying phase in different light windows.
%Can be used to optimize parameters for best fit or evaluate a given
%parameter set.

%2016-08-24, EL: test error analysis

%2016-08-23, EL: use to fit Vijayan et al. (2009, PNAS) data

%2015-10-21, EL
%Modified to simplify sine fit, removed sloping baseline
%Note: realized today that lsqcurvefit uses same algorithm as lsqnonlin,
%but the interface is simpler (returns yfit directly!). consider rewriting
%this.

%2015-09-20, EL
% Modified to return cell array params and fit ys
% pass in xdata, ydata as cell arrays; in function the params are recast as
% a matrix, and residuals returned as a vector (lsqnonlin requirement) but
% params and yfit are then repackaged again into cell arrays to be returned.
% Also modified to have the ability to pass in upper/lower bounds on period
% with ['periodUB', value] syntax (see parser code below).
% Also added the ability to evaluate a parameter set without doing any
% fits. Pass 'paramsToEval' argument, which must be a cell array or matrix
% of appropriate dimensions

%2015-03-20, EL
% This code fits multiple parts of a TopCount trace to indivdiual
% sinusoids sharing the same period but independent phase.
% datatofit = a cell array, each entry contains one "window" to be fit to a sinusoid
% timestofit = a cell array containing corresponding timepoints
% varargin = variable arguments possible, parser looks for period only
% fitresults = a matrix containing one set of fit params per window, each
% column contains [offset; lin. slope; period; amplitude; phase offset; chisqfit]

function [fitresults, fitdata, resnorm, jacobian] = ...
    fitSinusoid_Jac(timestofit, datatofit, varargin)
fitresults = [];
fitdata = [];
resnorm=[];
jacobian=[];

% set up input parser
p = inputParser;
defaultPeriod = 24;
addRequired(p,'timestofit',@iscell);
addRequired(p,'datatofit',@(x) (isnumeric(x) || iscell(x)));
%addOptional(p,'period',defaultPeriod,@(x) isnumeric(x));
addOptional(p,'period',defaultPeriod,@(x) x==x);
addOptional(p,'periodLB', [], @(x) (isnumeric(x) || isempty(x)));
addOptional(p,'periodUB', [], @(x) (isnumeric(x) || isempty(x)));
addOptional(p,'paramsToEval',[], @(x) (isnumeric(x) || iscell(x)));
addOptional(p,'paramsGuesses',[], @(x) (isnumeric(x) || iscell(x)));
parse(p,timestofit, datatofit, varargin{:});
periodToFit = p.Results.period;
periodLB = p.Results.periodLB;
periodUB = p.Results.periodUB;
paramsToEval = p.Results.paramsToEval;
paramsGuesses = p.Results.paramsGuesses;
% disp('paramsToEval:');
% disp(paramsToEval);

%if passed 'paramsToEval' function, then don't need to fit, just evaluate a
%parameter set
if ~isempty(paramsToEval)
    disp(['Evaluating parameters only, not fitting.']);
    NUMWINDOWS = numel(timestofit);
    paramsToEvalMat = [];
    chisq = []; 
    
    %if params are a cell array, convert to matrix
    %fitter requires a parameter matrix with each column = params for one
    %window
    if iscell(paramsToEval)
        assert(numel(paramsToEval) == numel(timestofit));
        for nw=1:NUMWINDOWS
            %disp(paramsToEval{nw});
            paramsToEvalMat(1:numel(paramsToEval{nw}),nw) = paramsToEval{nw}';
        end
    elseif isnumeric(paramsToEval) %if passing matrix, make sure it has right dimensions (NWINDOWS columns)
        paramsToEvalMat = paramsToEval; %must be a column vector
    else
        error('unrecognized parameter set type.');
    end
    
    %if passed datatofit, then want to return chisq
    if iscell(datatofit) && (numel(datatofit) == numel(timestofit))
        chisq = sinusoidSimple(paramsToEvalMat,timestofit,datatofit);
        chisq = sum(chisq.^2); %returns a vector, need to sum squares of residuals
        %disp(chisq);
    end
    
    %run sine fn, save ydata, chisq
    for nw=1:NUMWINDOWS
        %return row vectors
        if ~isempty(chisq)
            fitresults{nw} = [paramsToEvalMat(:,nw); chisq/(NUMWINDOWS*4+1)];
        else
            fitresults{nw} = [paramsToEvalMat(:,nw); -1];
        end
        fitdata{nw} = sinusoidSimple(fitresults{nw},timestofit(nw),0);
        fitresults{nw} = fitresults{nw}';
    end
    jacobian=[];
return;
end

%% if haven't returned, then must fit
assert(iscell(datatofit));

% options for fit
myoptions = optimoptions(@lsqnonlin,'maxfunevals',10^6); %,'InitDamping',100);
%myoptions = optimset('Largescale','on');
touseslope = 0;

% initialize guesses
NUMWINDOWS = numel(timestofit);

params0 = zeros(4, NUMWINDOWS);
lb = zeros(4, NUMWINDOWS);
ub = zeros(4, NUMWINDOWS);

% indices of fit parameters
offset = 1;
period = 2; 
amplitude = 3;
phase = 4;

% 1 - offset
params0(offset,:) = 0;  
lb(offset,:) = -100;
ub(offset,:) = 100;

% 2 - period, same fit for all windows
params0(period,1) = periodToFit; % use this entry for all fits
params0(period,2:end) = 0;
if isempty(periodLB)
    lb(period,:) = periodToFit - 2;
else
    lb(period, :) = periodLB;
end
if isempty(periodUB)
    ub(period,:) = periodToFit + 2;
else 
    ub(period,:) = periodUB;
end

% 3 - amplitude
params0(amplitude,:) = 100; 
lb(amplitude,:) = 0.0;
ub(amplitude,:) = 4*params0(amplitude,:);

% 4 - phase 
params0(phase,:) = 1.0*pi; 
lb(phase,:) = 0;
ub(phase,:) = 2*pi;

assignin('base','params0', params0);
assignin('base','lb', lb);
assignin('base','ub', ub);

disp(['Period passed to fitter: ' num2str(periodToFit)]);
disp(['period LB: ' num2str(lb(period,1)) ...
    ', period UB: ' num2str(ub(period,1))]);

[fit_params, resnorm, residual, exitflag, output ,lambda, jacobian] = ...
    lsqnonlin(@sinusoidSimple,params0,lb,ub,myoptions,timestofit,datatofit);

chisq=resnorm;
assignin('base', 'chisq', chisq);
 
% fitresults(1:NUMWINDOWS,1:5)=fit_params';
% fitresults(1:NUMWINDOWS,6)=chisq/(NUMWINDOWS*4+1); %reduced chisq. divide by num. fit. params.

%unpack params, get cell array fitys
%return row vectors
for i=1:NUMWINDOWS
    fitresults{i} = [fit_params(:,i); chisq/(NUMWINDOWS*4+1)];
    fitdata{i} = (sinusoidSimple(fitresults{i},timestofit(i),{zeros(size(timestofit(i)))}));
    fitresults{i} = fitresults{i}'; %must transpose after calling sinusoid_
end
assignin('base','jacobian',jacobian);
assignin('base', 'fitresults',fitresults);
assignin('base', 'fitdata',fitdata);
assignin('base', 'datatofit',datatofit);

end


