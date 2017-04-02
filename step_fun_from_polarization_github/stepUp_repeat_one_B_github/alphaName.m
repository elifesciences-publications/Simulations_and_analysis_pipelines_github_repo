%2017-02-13, EL
function [aName, newNum] = alphaName(numName)
%convert a number into base 26, alphabetic name, used by Excel in column
%names.

abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

assert(numName > 0);

newNum=[];
div = numName;
cnt=1;
while div > 0
    rem = mod(div, 26);
    div = (div-rem)/26;
    newNum(cnt) = rem;
    cnt=cnt+1;
end

%disp(newNum);

for n=1:numel(newNum)
    if newNum(n) < 1
        if n==numel(newNum)
            break;
        elseif newNum(n) == -1
            aName(n) = 'Y';
            newNum(n+1) = newNum(n+1)-1;
        elseif newNum(n) == 0
            aName(n) = 'Z';
            if n<numel(newNum)
                newNum(n+1) = newNum(n+1)-1;
            end
        end
    else
        aName(n) = abc(newNum(n));
    end
end

aName = fliplr(aName);
newNum = fliplr(newNum);

end

