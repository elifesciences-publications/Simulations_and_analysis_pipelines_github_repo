function [ OUTPAR_CA ] = convVecParToCA(INPAR_MAT)
%turn row of 31 parameters into a cell array. 31 params are arranged in 10
%sets of [offset amplitude phase]x10 followed by 31st entry which is period
%shared by all 10 sets. Return 10x1 cell array with each cell having
%{offset, period, amplitude, phase}x10, ordered in the same way as the 
%parameters in the row vector.
%if INPAR is a matrix, treat every row in this way. output a cell array
%where each entry is also a cell array with 10x4 cells.

[numrow numcol] = size(INPAR_MAT);

%each row must have 31 params
assert(numcol == 31);

for n=1:numrow
    pars = INPAR_MAT(n,:);
    
    par_ca = [];
    for r=1:10
        %[1 31 2 3], [4 31 5 6], [7 31 8 9], etc.
        par_ca{r} = pars([3*r-2 31 3*r-1 3*r]);
    end
    
    OUTPAR_CA{n} = par_ca;
end

end

