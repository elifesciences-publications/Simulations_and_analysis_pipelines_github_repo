%EL, 2016-10-26. Modify to weigh residuals based on errors on observations,
%i.e., compute chi_i = (yfit_i-y_i)/yerr_i at every point.
%EL, 2016-04-28. Fit multiple datasets to multiple functions that share
%function parameters.
%Inputs: params = [coefs]; DATA and FITFUNS are structures that have the
%same number of objects. For every dataset, DATA contains a multi-column
%matrix where dependent variables are in last column; FITFUNS contains 
%handles to functions to which you're fitting the data.
%Outputs: [result] is a vector of residuals. Computed point-by-point by
%default.

function [result] = globalMultiLinFitFun(params,DATA,FITFUNS)

disp(mfilename('fullpath'));

TOTEST = 0;
if TOTEST == 1
    data.a = [1 2 3; 4 5 6; 7 8 9];
    data.b = [1 3; 4 0; 12 10; 34 20];
    
    fitfun.a = @(x) x(:,end-1) + 1;
    fitfun.b = @(x) x(:,end-1) - 1;
    
    DATA = data;
    FITFUNS = fitfun;
end

    setnames = fieldnames(DATA);
    funnames = fieldnames(FITFUNS);
    assert(numel(setnames) == numel(funnames));
    
    numsets = numel(setnames);
        
    result = [];
    
    for n=1:numsets
         thisData = DATA.(setnames{n});
         thisFun = FITFUNS.(funnames{n});
%          disp(thisData);
%          disp(thisFun);
         disp([setnames{n} ' *** ' funnames{n}]);
       
         thisFit = thisFun(params, thisData);
         
         %assume ydata in col 3, yerr in col 4; ystd in col. 5
         thisRes = (thisData(:,3) - thisFit)./thisData(:,5); %divided by thisData(:,5) in fig files
         
         %should be adding residual inside the loop
         %try scaling residuals by no. pts: (res/numel(thisFit))
         result = [result; thisRes];
    end
       
end

