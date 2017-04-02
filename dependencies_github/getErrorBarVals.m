function [xout, yavg, yerr, ystd, ynum ] = getErrorBarVals(x, y)
%Pass in a list of xs and corresponding ys. For every unique x, compute avg
%y, yerr, ystd and numpts. Useful for subsequent errorbar() call. Excludes
%any NaN values in y by calling nanmean, nanstd. For NaN values in x, 
%simply exludes all [x y] pairs where x is NaN.
%IN: [x], [y]. OUT: [xout, yavg, yerr, ystd, ynum],

xout=[];
yavg=[];
yerr=[];
ynum=[];

assert(numel(x) == numel(y));
xout = unique(x);
xout = xout(~isnan(xout));
%disp(xout);

for i=1:numel(xout)
    thisx = xout(i);
    xind = (x == thisx);
    thisy = y(xind);
    yavg(i) = nanmean(thisy);
    numnans = sum(~isnan(thisy));
    yerr(i) = nanstd(thisy)/sqrt(numnans);
    ystd(i) = nanstd(thisy);
    ynum(i) = numnans;
end

end
