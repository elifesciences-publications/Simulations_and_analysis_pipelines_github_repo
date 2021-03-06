%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook:
%    /Users/E/Documents/Advisers/Rust/Clocks/2015-11-23_stepUp-stepDown/2015-12_08_sU_quants/2015-12-08_sU_quants.xlsx
%    Worksheet: Sheet1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2015/12/08 18:17:34

%% Import the data
function [Rxn, TubeTime, ActTime, pSpT, pT, pS, pU, pNet] = ...
    loadQuants(FILENAME)
[~, ~, raw] = xlsread(FILENAME, 'data');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
raw = raw(:,[2,3,4,5,6,7,8]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
Rxn = cellVectors(:,1);
TubeTime = data(:,1);
ActTime = data(:,2);
pSpT = data(:,3);
pT = data(:,4);
pS = data(:,5);
pU = data(:,6);
pNet = data(:,7);

%% Clear temporary variables
clearvars data raw cellVectors R;