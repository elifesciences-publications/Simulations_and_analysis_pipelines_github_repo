function [lincoef, coef_err] = fitToLineYErr(fitxs, fitys, yerr)
%fit portions of curve to lines given [x, y, yerr] vectors as inputs
%return ([slope yin], [slope-err yint-err adj_rsq])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %NR assumes fit is of the form:
    %    y(x; a,b) = a + b*x
    %from numerical recipes 
    %   (15.2; eq. 15.2.4-9)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S = sum(1./(yerr.^2));
    Sx = sum(fitxs./(yerr.^2));
    Sy = sum(fitys./(yerr.^2));
    Sxx = sum((fitxs.^2)./(yerr.^2));
    Sxy = sum((fitxs.*fitys)./(yerr.^2));
    
    D = S*Sxx - (Sx)^2;
    
    a = (Sxx*Sy - Sx*Sxy)/D;
    b = (S*Sxy-Sx*Sy)/D;
    
    a_sigma = sqrt(Sxx/D);
    b_sigma = sqrt(S/D);

    lincoef = [b a];
    coef_err = [b_sigma a_sigma ];
    
    %compute adj rsq using old formula (Matlab documentation)
    linfit =  polyval(lincoef,fitxs);
    resid = fitys - linfit;
    SSresid = sum(resid.^2);
    SStotal = (length(fitys)-1)*var(fitys);
    rsq_adj = 1 - SSresid/SStotal * ...
        (length(fitys)-1)/(length(fitys)-length(lincoef));
    
    coef_err = [coef_err rsq_adj];
end