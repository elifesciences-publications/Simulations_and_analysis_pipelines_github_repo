function [ OUTPAR_MAT ] = convCAParToVec(INPAR_CA)
%convert parameter sets in cell array form into a matrix.
%   Each cell of a cell array is one parameter set for 10 reactions. Each
%   parameter set itself is a cell array of ten 1x4 vectors, each contaning
%   [offset period amplitude phase] for one curve. Convert them into a
%   matrix, where each row is one parameter set, ordered the same way as
%   they are ordered in the cell array, but having the following pattern
%   10x[offset amplitude phase] then [period], since period is fit globally.

%if you only have 1 parameter set, enclose it in another layer of {}: {{[],[],...[]}}
numrow = numel(INPAR_CA);

for n=1:numrow
  assert(numel(INPAR_CA{n}) == 10);
  parset = INPAR_CA{n}; %should be a cell array with 10 entries
    
  parvec = [];
  for r=1:10
    parvec([3*r-2 3*r-1 3*r]) = parset{r}([1 3 4]);
    if r==1
        parvec(31) = parset{r}(2);
%         disp('set parvec31 = '); 
%         disp(parvec(31));
    else
%         disp('new parvec31 attempt=');
%         disp(parset{r}(2));
        %period should be the same in all rxns
        assert(parset{r}(2) == parvec(31)); 
    end
  end
  
  OUTPAR_MAT(n,:) = parvec;
  
end

end

