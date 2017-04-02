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

function [result] = sinusoidSimple(params,X,Y)
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
        period = params(2,1); % always use first entry for period
        amplitude = params(3,i);
        phase = params(4,i);
        
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
end

