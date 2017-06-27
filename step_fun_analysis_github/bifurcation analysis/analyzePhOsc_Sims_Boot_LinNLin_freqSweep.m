%2017-03-28, EL: Fig. 3---figSup2 
%process bifurcation simulation results and make a bifucation diagram

clear all;
close all;
clc;
cd(['.']);
%% load simulation results

SIMULATION_FILE = load('2017-06-05_14.01.10_freqSweep_stepFunMixMatch_NonLin_1.mat');
SIM = SIMULATION_FILE.VAR;
dawnphase = SIM.dawnphase;
duskphase = SIM.duskphase;
duskphaseshift = SIM.duskphaseshift;
dawnphaseshift = SIM.dawnphaseshift;
peakT = SIM.peakT;

%% simulation parameters
TDRIVE=[6:0.0025:48];
numcyc = [1:100*100]; %should be a vector, not just max value!

DUTYFRAC=0.5; %always balanced duty cycle
TEND = (max(numcyc)+0.1)*TDRIVE;

%starting phases for each dutyfrac
startphase = zeros(size(TDRIVE)); %[0];

%% collect variables for plotting
plotphase = dawnphase; 
plotphase = mod(plotphase,1);

%% plot all visited dawn phases as a fun of driving period
plotphase = mod(plotphase,1);

fSimpleBifurc = figure();

for fr=1:numel(TDRIVE)
    longphases = plotphase(fr,50:end);
    driveper = ones(size(longphases))*TDRIVE(fr);
    plot(driveper,longphases, 'k.', 'markersize', 1);
    hold on;
    plot(driveper,1+longphases, 'k.', 'markersize', 1);
end

grid off;
set(gca, 'xlim', [0 50], 'xtick', [0:4:48], 'ylim', [0 2], 'ytick', [-2:0.2:2]);
xlabel('driving period (hours)');
ylabel('phase at end of nighttime (rad/2\pi)');

set(gcf,'units','inches','position',[0 0 10 5]);
%%
TOEXP_BIFURC = 1;
if TOEXP_BIFURC == 1
    FILENAME = [getDate('yyyy-mm-dd') ...
        '_bifurcSimple_Lin-th0-' num2str(startphase(1)) '_' ...
        getDate('HH.MM.SS')];
   
    %looks like 'r1000' flag generates a png with enough pixel density that
    %you can zoom in a few times without the image looking terrible. fonts
    %also look ok.
    export_fig(FILENAME,'-cmyk', '-r1000', '-nocrop', '-jpg', fSimpleBifurc);
end

break;
