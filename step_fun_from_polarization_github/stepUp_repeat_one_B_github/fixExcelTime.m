%2017-02-13, EL
function [fullDateVec, fullDateStr, fullDateNum] = fixExcelTime(vecDate, txtTime)
%Convert between Excel and Matlab time formats. Pass in date as a datevec, 
%time as txt. This input assumes that Tecan iControl software records the
%date as a 'DD-MMM-YYYY' and time as 'HH:MM:SS AMPM'.

cellDate = strsplit(datestr(vecDate),'-');
d = str2num(cellDate{1}) - 1;
m = cellDate{2};
y = 1900 + str2num(cellDate{3});

if numel(d) < 2
    strDate = ['0' num2str(d) '-' m '-' num2str(y)];
else
    strDate = ['0' num2str(d) '-' m '-' num2str(y)];
end

dVec = datevec(strDate);
tVec = datevec(txtTime);

fullDateVec = [dVec(1:3) tVec(4:6)];
fullDateStr = datestr(fullDateVec);
fullDateNum = datenum(fullDateVec);


end

