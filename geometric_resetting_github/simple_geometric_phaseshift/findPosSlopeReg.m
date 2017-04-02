function [ YOUT, IXOUT ] = findPosSlopeReg(YIN)
%Given a vector YIN, return the longest continuous portion of YIN that has
%a positive slope throughout. Return vector YOUT of values in addition to
%their indices in the original YIN vector.

%tested with: test = [-1 0 2 3 0 0 1 2 3 2 2 3 4 5]; 

y = YIN;
dy = diff(y,1);

%find where 1st derivative changes sign
pos_dy = dy > 0;

%go thru list
startpt = find(pos_dy,1,'first');
endpt = startpt+1;
n = 1; %number of positive runs
runs = [];

while endpt <= numel(y)
    if y(endpt) > y(endpt-1)
        %always save last run
        if endpt == numel(y)
            runs(n,1:2) = [startpt endpt];
            endpt = endpt+1;
        else
            endpt = endpt+1;
        end
    else
        endpt=endpt-1;
        runs(n,1:2) = [startpt endpt];
        
        pos_dy = dy(endpt+1:end) > 0;
        startpt = endpt+find(pos_dy,1,'first');
        endpt = startpt+1;
        n = n+1;
    end
end

%if don't have any runs > 1
if isempty(runs)
    runs(1,:) = [1 1];
end

runs(:,3) = runs(:,2) - runs(:,1);
[~,maxrun] = max(runs(:,3));
IXOUT = runs(maxrun,1):runs(maxrun,2);
YOUT = [y(IXOUT)];

end

