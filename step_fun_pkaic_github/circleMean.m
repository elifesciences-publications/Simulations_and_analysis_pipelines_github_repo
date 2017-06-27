function [MEANOUT, STDEVOUT, OUTBREAKPT] = circleMean(VEC, PER, INBREAKPT)
%Take the mean and stdev of data on the circle. VEC = vector of data on the
%circle, PER = period of the data or circumference of the circle, INBREAKPT =
%phase (0-PER) to be used as phase0 (e.g., midnight for a clock). Pass
%INBREAKPT = [] to determine the breakpt automatically. 
%Return MEANOUT, STDEVOUT of data on the circle as well as the 
%computed breakpoint as OUTBREAKPT. Outputs in same coordinates as VEC (not
%in phase coordinates).

%get rid of nan's for simplicity
vec = VEC(~isnan(VEC));

if numel(vec) < 1
    warning('input vec is all NaN!');
end

%convert data to polar coordinates: x=r*cos(theta), y=r*sin(theta)
r = PER/(2*pi);
theta = (2*pi)*mod(vec,PER)/PER;
x = r.*cos(theta);
y = r.*sin(theta);

%find mean phase of datapts
if ~isempty(INBREAKPT)
    thetabar=INBREAKPT;
else
    xbar = mean(x);
    ybar = mean(y);
    thetabar = atan2(ybar, xbar); %find breakpt
end

%set breakpoint, arrange data on [-pi, pi] relative to breakpt
vecprime = mod(theta - thetabar, 2*pi);
vecprime(vecprime > pi) = vecprime(vecprime > pi) - 2*pi;
vecprime(vecprime <= -pi) = vecprime(vecprime <= -pi) + 2*pi;

%compute mean using breakpoint as reference
MEAN = nanmean(vecprime) + thetabar;
STD = nanstd(vecprime);
BREAKPT = thetabar;

%return values in same coordinates as inputs
MEANOUT = MEAN*PER/(2*pi);
STDEVOUT = STD*PER/(2*pi);
OUTBREAKPT = BREAKPT*PER/(2*pi);

if isnan(MEANOUT) || isnan(STDEVOUT)
    warning('Problem in circle mean!');
    disp(['MEAN: ' num2str(MEAN) ', STD: ' num2str(STD) ...
          ', break: ' num2str(BREAKPT) ', PER: ' num2str(PER) ...
          ', thetabar: ' num2str(thetabar)]);
end

%check by plotting
TOPLOT = 0;
if TOPLOT == 1
   fTest = figure();
   t=[0:0.01:2*pi];
   
   % plot circle of radius 1
   plot(r*sin(t), r*cos(t), 'k-');
   hold on;
   
   %plot input pts
   plot(x,y,'bo');
   
   %plot avg
   plot(xbar,ybar,'rs');
   plot(r*cos(MEAN), r*sin(MEAN), 'go');
   plot(r*cos(MEAN-STD), r*sin(MEAN-STD),'gx');
   plot(r*cos(MEAN+STD), r*sin(MEAN+STD),'gx');
   axis('square');
end

end

