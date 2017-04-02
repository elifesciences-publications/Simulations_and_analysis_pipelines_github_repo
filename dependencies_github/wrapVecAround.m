function [WRAPMAT] = wrapVecAround(MAT, VAL, PERIOD, SIGN)
%Pass matrix MAT (or vector) and value VAL around which you want to wrap
%all values. If SIGN = 'gt' or 'greaterthan', mat(mat > val) = mat(mat > val) - period;
%if SIGN = 'lt' or 'lessthan', mat(mat < val) = mat(mat < val) + period.

if strcmp(SIGN, 'gt') || strcmp(SIGN,'greaterthan')
    MAT(MAT > VAL) = MAT(MAT > VAL) - PERIOD;
elseif strcmp(SIGN,'lt') || strcmp(SIGN, 'lessthan')
    MAT(MAT < VAL) = MAT(MAT < VAL) + PERIOD;
else
    warning(['did not pass sign; no wrapping done']);
end

WRAPMAT = MAT;

end

