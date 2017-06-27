%2017-03-30, EL: generate Fig. 4-figSup2(A-B).
%Plot boostrapped distributions of L and D funs as error bars with x and y
%errors.
%
% Enter the file name containing bootstrapped L and D step functions in
% INFILE_FUN variable.
% 
% dependencies: errorbarxy.m from FileExchange to generate errorbar plots,
%               export_fig.m from FileExchange for figure export.
%

close all;
clear all;
%% input parameters
TOEXP_STEPFUNS = 0; %export fig?

INFILE_FUN = 'saved_data/2017-06-05_stepFun_JuneFit_updFitPipe_wideLinFit1000_15.49.50.mat'; %
load(INFILE_FUN); %stepFuns and fit periods

%helper fun
saveMyFig = (@(fig,figName) ...
    export_fig([getDate('yyyy-mm-dd') '_' figName '_' getDate() '.pdf'],...
    '-cmyk', '-painters','-pdf',fig));

%% deal with stepFuns: get x, y errorbars 

for np=1:numel(up.phase(1,:))
    [mean_upPhase(np), std_upPhase(np), up_break(np)] = ...
        circleMean(up.phase(:,np),4*2*pi,[]);
    [mean_upPhaseShift(np), std_upPhaseShift(np), upShift_break(np)] = ...
        circleMean(up.phaseShift(:,np),4*2*pi,[]);
end

for np=1:numel(down.phase(1,:))
    [mean_downPhase(np), std_downPhase(np), down_break(np)] = ...
        circleMean(down.phase(:,np),4*2*pi,[]);
    [mean_downPhaseShift(np), std_downPhaseShift(np), downShift_break(np)] = ...
        circleMean(down.phaseShift(:,np),4*2*pi,[]);
end
%% wrap stepDown around 2*pi
mean_upPhase = mean_upPhase + pi/2; %add pi/2 to count relative to the min of %P KaiC trajectory
% mean_upPhase = wrapVecAround(mean_upPhase,4*pi,4*pi,'gt');
mean_downPhase = mean_downPhase + pi/2; % add pi/2 to count relative to the min of %P KaiC trajectory
% mean_downPhase = wrapVecAround(mean_downPhase,4*pi,4*pi,'gt');

%sort stepfuns
[~,downIx] = sort(mean_downPhase);
[~,upIx] = sort(mean_upPhase);

mean_upPhase = mean_upPhase(upIx);
mean_upPhaseShift = mean_upPhaseShift(upIx);
std_upPhase = std_upPhase(upIx);
std_upPhaseShift = std_upPhaseShift(upIx);

mean_downPhase = mean_downPhase(downIx);
mean_downPhaseShift = mean_downPhaseShift(downIx);
std_downPhase = std_downPhase(downIx);
std_downPhaseShift = std_downPhaseShift(downIx);

%% plot stepup, use errorbarxy(x,y,dx,dy)
z = (24)/(2*pi);

fUpNoisy = figure();
[hdx_up,hdy_up] = errorbarxy(z*mean_upPhase,z*mean_upPhaseShift, ...
    z*std_upPhase, z*std_upPhaseShift,'linewidth',1,'color','b');
hold on;
plot(z*mean_upPhase, z*mean_upPhaseShift,'sb-',...
    'markersize',4, 'markerfacecolor','b'); %,...
    %'markeredgecolor','none','linewidth',1);
grid off;
% uistack(hdx_up,'top');
% uistack(hdy_up,'top');
box on;
set(gca,'xlim',z*[0*pi 4*pi],'xtick',z*[-2*pi:0.5*pi:4*pi],...
    'ylim',z*[-1*pi 1.3*pi],'ytick',z*[-pi:(pi/3):pi]);
xlabel('step time (CT hours)');
ylabel('phase shift (CT hours)');

fDownNoisy=figure();
[hdx_down,hdy_down] = errorbarxy(z*mean_downPhase, z*mean_downPhaseShift, ...
    z*std_downPhase, z*std_downPhaseShift,'linewidth',1,'color','r');
hold on;
plot(z*mean_downPhase, z*mean_downPhaseShift,'sr-',...
    'markersize',4, 'markerfacecolor','r'); %,...
    %'markeredgecolor','none','linewidth',1);
grid off;
% uistack(hdx_up,'top');
% uistack(hdy_up,'top');
set(gca,'xlim',z*[0*pi 4*pi],'xtick',z*[-2*pi:0.5*pi:4*pi],...
    'ylim',z*[-pi 1.5*pi],'ytick',z*[-pi:(pi/3):pi]);
box on;
xlabel('step time (CT hours)');
ylabel('phase shift (CT hours)');

set(fDownNoisy,'units','inches','position',[3 3 3.75 3]);
set(fUpNoisy,'units','inches','position',[3 3 3.75 3]);

%% add in vitro data
%label points based on when dusk/dawn hit in ATP/ADP driving (Fig. 2C)

%this data is in [0,1] phase units. aligned such that ph0 = min of KaiC
%phosphorylation
inVitroDawnPhase =[  0.1420    0.2507    0.1820    0.1450    0.0958    0.0619   -0.0131];
inVitroDuskPhase =[  1.1441    0.5002    0.5146    0.5608    0.5969    0.6465    0.7385];
inVitroPP =[24 6 8 10 12 14 18];

colors = hot(numel(inVitroDawnPhase));
dusk = mod(2*pi*z*inVitroDuskPhase,2*pi*z);
dawn = mod(2*pi*z*inVitroDawnPhase,2*pi*z);
dawn = wrapVecAround(dawn,4*pi*z,2*pi*z,'gt');
dawn = wrapVecAround(dawn,pi*z,2*pi*z,'lt');

duskStep = interp1(z*mean_downPhase, z*mean_downPhaseShift,dusk);
dawnStep = interp1(z*mean_upPhase, z*mean_upPhaseShift,dawn);
labels = {'24','6:18','8:16','10:14','12:12','14:10','18:16'};

%%
% for i=2:numel(dusk)
%     figure(fDownNoisy);
%     plot(dusk(i),duskStep(i),'o','markerfacecolor',colors(i-1,:),...
%         'markersize',4,'markeredgecolor','k','linewidth',0.5);
%     figure(fUpNoisy);
%     plot(dawn(i),dawnStep(i),'o','markerfacecolor',colors(i-1,:),...
%         'markersize',4,'markeredgecolor','k','linewidth',0.5)
% end

for i=2:numel(dusk)
    figure(fDownNoisy);
    z=24/(2*pi); %this z must be in the same units as in the plots above
    
    %try arrows
    arrow([dusk(i) -z*0.7*pi],[dusk(i) -z*0.45*pi],'length',8,...
        'facecolor',colors(i,:),'edgecolor','k','tipangle',20);
    hold on;
    
    figure(fUpNoisy);
    %try arrows
    arrow([dawn(i) -z*0.7*pi],[dawn(i) -z*0.45*pi],'length',8,...
        'facecolor',colors(i,:),'edgecolor','k','tipangle',20);
    hold on;
end

%% export data into CSV files
dsUpKaiCP = z*[mean_upPhase' std_upPhase' mean_upPhaseShift' std_upPhaseShift'];
dsUpKaiCP_names = {'step_up_time','step_up_time_err','step_up_phaseshift','step_up_phaseshift_err'};
dsUpKaiCP = mat2dataset(dsUpKaiCP);
dsUpKaiCP.Properties.VarNames=dsUpKaiCP_names;

dsDownKaiCP = z*[mean_downPhase' std_downPhase' mean_downPhaseShift' std_downPhaseShift'];
dsDownKaiCP_names = {'step_down_time','step_down_time_err','step_down_phaseshift','step_down_phaseshift_err'};
dsDownKaiCP = mat2dataset(dsDownKaiCP);
dsDownKaiCP.Properties.VarNames=dsDownKaiCP_names;

dsDusk = mat2dataset([inVitroPP(2:end)' dusk(2:end)' duskStep(2:end)']);
dsDusk_names = {'day_length', 'dusk_time_CT', 'dusk_phaseshift_CT'};
dsDusk.Properties.VarNames = dsDusk_names;

dsDawn = mat2dataset([inVitroPP(2:end)' dawn(2:end)' dawnStep(2:end)']);
dsDawn_names = {'day_length', 'dawn_time_CT', 'dawn_phaseshift_CT'};
dsDawn.Properties.VarNames = dsDawn_names;

% export(dsUpKaiCP,'FILE','L_function.csv','Delimiter',',');
% export(dsDownKaiCP,'FILE','D_function.csv','Delimiter',',');
% export(dsDawn,'FILE','dawns.csv','Delimiter',',');
% export(dsDusk,'FILE','dusks.csv','Delimiter',',');

%% export temporarily to .mat files
% save([getDate('yyyy-mm-dd') '_PKaiCfuns_' getDate('HH.MM.SS') '.mat'],...
%     'dsUpKaiCP','dsDownKaiCP');

%% export step fun figures
if TOEXP_STEPFUNS == 1
    export_fig([getDate('yyyy-mm-dd') '_NPboot_stepFun_ErrBars_up-' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',fUpNoisy);
    export_fig([getDate('yyyy-mm-dd') '_NPboot_stepFun_ErrBars_down-' getDate() '.pdf'],...
        '-cmyk','-painters','-pdf',fDownNoisy);
end

