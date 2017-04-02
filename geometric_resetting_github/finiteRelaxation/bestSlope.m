function [slope, rsq, breakpt, goodfit, best_err] = ...
    bestSlope(xin, yin, per, posneg,varargin)
%Given a vector [x1 y1; x2 y2; ...] where y's are periodic with period per,
%find the best fit slope to the points, taking periodicity into account. In
%other words, make sure all of the points are within one period of each
%other and that adding or subtracting a period from any point isn't going
%to change the value of the slope.
%Return best-fit slope, rsq and breakpoint for yin; goodfit = 0 or 1
%depending on rsq.

TOTEST = 0;
if TOTEST == 1
    %test vals
    posneg = -1;
    xin = [6.1:3:18.1];
    per = 24;
    y = -0.7*xin + 10;
    yoffset = [per -per +2*per 0 0 0 per];
    ynoise = 2*randn(size(y));
    yin = y + ynoise;  
end

%parse input to determine if you want to use rsq or sqerr to judge fit
%quality. SQERR by default.
USE_SQERR = 1;
USE_RSQ = 0;
if ~isempty(varargin)
    assert(numel(varargin) == 1 && ischar(varargin{1}));
    switch varargin{1}
        case 'rsq'
            USE_SQERR = 0;
            USE_RSQ = 1;
        case 'sqerr'
            disp('using squerr');
            USE_SQERR = 1;
            USE_RSQ = 0;
    end
end

%sort inputs
[~,ix] = sort(xin);
xin=xin(ix);
yin=yin(ix);

%make sure all ys are within one period of each other
% while max(yin) - min(yin) > per
%     [maxy,maxy_ix] = max(yin);
%     yin(maxy_ix) = yin(maxy_ix) - per;
% end
yorig = yin;
yin = mod(yin,per);


%do linear fits using every yin value as a breakpoint
%try to enforce linearity by placing points to the left < X, to the right >
%X, by taking advantage of periodicity

%if ascending
for n=1:numel(yin)
    ywrap = yin;
    ywrap((n+1):end) = wrapVecAround(ywrap((n+1):end),yin(n),per,'lt');
    ywrap(1:n) = wrapVecAround(ywrap(1:n),yin(n),per,'gt');
    [fitval(n,:) fitrsq(n)] = fitToLine(xin,ywrap); 
    squerr(n) = sum((polyval(fitval(n,:),xin) - ywrap).^2);
end

N = numel(yin);
%if descending
for n=1:numel(yin)
    ywrap = yin;
    ywrap((n+1):end) = wrapVecAround(ywrap((n+1):end),yin(n),per,'gt');
    ywrap(1:n) = wrapVecAround(ywrap(1:n),yin(n),per,'lt');
    [fitval(N+n,:) fitrsq(N+n)] = fitToLine(xin,ywrap); 
    squerr(N+n) = sum((polyval(fitval(n,:),xin) - ywrap).^2);
end

%try fitting without scrambling the order
[fitval(2*N+1,:) fitrsq(2*N+1)] = fitToLine(xin,yin); 
squerr(2*N+1) = sum((polyval(fitval(2*N+1,:),xin) - yin).^2);

%%%%%%%%%% PROBLEM !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MUST DO THIS SEPARATELY FOR DOWNWARD & UPWARD SLOPING POINTS
%%% MESSY. TRY CIRCLEMEAN STRATEGY INSTEAD?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%decide if slope has to be positive, neg. or either and whether to use rsq
%or sqerr to judge fit quality
if USE_RSQ == 1
    if posneg == 1
        [~, bestix] = max(fitrsq(1:N));
    elseif posneg == -1
        [~, bestix] =  max(fitrsq((N+1):end));
        bestix = N + bestix;
    else
        [~, bestix] = max(fitrsq);
    end
elseif USE_SQERR == 1
    if posneg == 1
        [~, bestix] = min(squerr(1:N));
    elseif posneg == -1
        [~, bestix] =  min(squerr((N+1):end));
        bestix = N + bestix;
    else
        [~, bestix] = min(squerr);
    end
end

slope = fitval(bestix,1);
rsq = fitrsq(bestix);
best_err = squerr(bestix);

%disp(bestix);
if bestix == 2*N+1
    breakpt = nan;
elseif bestix > N && bestix < 2*N + 1
    breakpt = yin(bestix-N);
else
    breakpt = yin(bestix);
end

%get pts in order used for fitting
y2fit = yin;

if bestix <= N
    y2fit((bestix+1):end) = ...
        wrapVecAround(yin((bestix+1):end),yin(bestix),per, 'lt');
    y2fit(1:bestix) = wrapVecAround(yin(1:bestix),yin(bestix),per,'gt');
elseif bestix < 2*N+1
    bx = bestix-N;
    y2fit((bx+1):end) = ...
        wrapVecAround(yin((bx+1):end),yin(bx),per, 'gt');
    y2fit(1:bx) = wrapVecAround(yin(1:bx),yin(bx),per,'lt');
end

%check fit quality
avgdev = mean(abs(polyval(fitval(bestix,:),xin) - y2fit));
if avgdev > 0.1*(max(y2fit)-min(y2fit)) && avgdev > 0.5 %best_err > 0.5 %place relative and absolute bounds
    goodfit = 0;
else
    goodfit = 1;
end

%check by plotting
TOPLOT = 0;
if TOPLOT == 1
    f=figure();
    %plot raw xin, yin
    plot(xin,yin, 'bo', 'linewidth',1);
    hold on;
    plot(xin,yorig,'kx','linewidth',1);
    
    %plot yvals used to fit
    plot(xin,y2fit,'ro');
      
    %plot fit
    plot(xin, polyval(fitval(bestix,:),xin),'k--');
    
    %set y axis limits
    minY = min([min(yin) min(yorig) min(polyval(fitval(bestix,:),xin))]);
    maxY = max([max(yin) max(yorig) max(polyval(fitval(bestix,:),xin))]);
    if maxY < minY + 24
        maxY = minY + 24;
    end
    set(gca, 'ylim', [minY maxY]);
    
    title(['y=' num2str(slope,'%2.2f') 'x + ' num2str(fitval(bestix,2)) ...
        '. R^2 = ' num2str(rsq,'%2.2f') ...
        '. sqErr = ' num2str(best_err,'%2.2f') ...
        '. goodfit = ' num2str(goodfit)]);
end


end

function [lincoef, rsq_adj] = fitToLine(fitxs, fitys)
%fit portions of curve to lines given [x, y] vectors as inputs
%return [slope y-int adj_rsq]
    lincoef = polyfit(fitxs,fitys,1);
    linfit =  polyval(lincoef,fitxs);
    resid = fitys - linfit;
    SSresid = sum(resid.^2);
    SStotal = (length(fitys)-1)*var(fitys);
    rsq_adj = 1 - SSresid/SStotal; %
                                %* (length(fitys)-1)/(length(fitys)-length(lincoef));
end




