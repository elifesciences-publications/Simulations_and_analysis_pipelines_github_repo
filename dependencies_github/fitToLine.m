function [lincoef, rsq_adj] = fitToLine(fitxs, fitys)
%fit portions of curve to lines given [x, y] vectors as inputs
%return [slope y-int adj_rsq]
    lincoef = polyfit(fitxs,fitys,1);
    linfit =  polyval(lincoef,fitxs);
    resid = fitys - linfit;
    SSresid = sum(resid.^2);
    SStotal = (length(fitys)-1)*var(fitys);
    rsq_adj = 1 - SSresid/SStotal * (length(fitys)-1)/(length(fitys)-length(lincoef));
end