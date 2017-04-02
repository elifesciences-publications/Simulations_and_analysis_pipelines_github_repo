%EL, 2016-09-21: use fixed period (=24 hrs).
%EL, 2016-09-21: use with fitEz.m
%EL, 2015-10-21. Modified to y=offset+ampl.*sin(2*pi*x/T - ph), remove slope.
%EL, 2015-09-20. Modified to include cell array outputs. 
%EL, 2015-03-20. Compute residual between a sinusoid fit and data for
% multiple disjointed intervals. All intervals share the same period, but
% phase, amplitude, offset change for diff. intervals,
% 
% X = cell array, each entry contains times for one lights-on window
% Y = cell array, each entry contains data to fit to sine for one lights-on
% window
% params = matrix, Nth column is a set of parameters to fit the Nth
% lights-on window. Note that all windows share the same period.
% order of params: [offset; period; amplitude; phase].

function [result, Jac] = ezSine(params,X,Y)
    %disp(['in multi-window sinusoid fn']);
    assignin('base','X',X);
    assignin('base','Y',Y);
    assignin('base','params_insine',params);    

    if iscell(X) && iscell(Y)
        NUMWINDOWS = numel(X);
    elseif isnumeric(X) && isnumeric(Y)
        NUMWINDOWS = 1;
    elseif Y==0
        NUMWINDOWS = numel(X);
        %warning('Y is 0!')
    else
        %warning('X and Y of different types!');
    end
    result = [];
  
    for i=1:NUMWINDOWS
        offset = params(1,i);
        period = 24; %params(2,1); % always use first entry for period
        amplitude = params(2,i);
        phase = params(3,i);
        
        %get xs, ys, accomodate empty ys
        xs = X{i};
        if iscell(Y)
            ys = Y{i};
        else
            ys = zeros(size(xs));
        end
        
        residual = offset + ...
            amplitude.*sin(2*pi.*(xs./period)- phase) - ys;
        
        % result needs to have more points than the no. of parameters
        % you're fitting; so keep residuals for all datapoints together
        % in a row vector
        %result = [result residual'];
        result = [result residual];     
    end
    assignin('base', 'result_insine', result);
    
    
    %compute jacobian if two output args
    if nargout > 1   
        Jac = [];   % Jacobian of the function evaluated at x
        %First column J is derivative of [result] with respect to params(1);
        %Second column of J is derivative of [result] wrt params(2), etc.
        
        dOffset = ones(size(result));
        dA = offset + sin(2*pi.*(xs./period)- phase);
%         dT = amplitude.*cos(2*pi.*(xs./period)- phase).*...
%             (-2*pi.*xs/(period.^2));
        dPhi = amplitude.*cos(2*pi.*(xs./period)- phase) -1;
        
        Jac = [dOffset dA dPhi];
        disp(Jac);
    end

end

%% From lsqnonlin documentation:
%If fun returns a vector (matrix) of m components and x has length n, where
%n is the length of x0, the Jacobian J is an m-by-n matrix where J(i,j) is
%the partial derivative of F(i) with respect to x(j). (The Jacobian J is
%the transpose of the gradient of F.)

%dRes/dOffset = 1;
%dRes/dA = offset + sin(2pix/T-ph)
%dRes/dT = Acos(2pi*x/T - ph)*(-2*pi*x/T^2)
%dRes/dPhi = Acos(2pi*x/T-ph)-1
