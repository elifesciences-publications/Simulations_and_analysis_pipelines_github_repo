%2017-03-28, EL: make Figure 3C-D.
%plot all experimentally measured L and D functions and mark with arrows where
%dawn and dusk occur.

%close all; 
clear all; clc;

INDIR = ['.'];
cd(INDIR);

%% parameters
TOEXP = 0; %export figure? (1=yes)
TOSAVEMERGE = 0; %save merged L and D datasets? (1=yes)

%% load all step funs from Feb - Mar measurements
%merge all step-response function measurements
STEPUP_ONE_DIR = '../stepUp_repeat_one_B_github/saved_data';
STEPUP_ONE_FILE = '2017-04-01_stepFun_Feb_plateReader_Apr1analyze_1_2017-04-01_11.29.55.mat';
step_up_one = load([STEPUP_ONE_DIR '/' STEPUP_ONE_FILE]);
step_up_one.up.name = 'step_up_one';
step_up_one.down.name = 'step_up_one';

STEPDOWN_ONE_DIR = '../stepDown_repeat_one_github/saved_data';
STEPDOWN_ONE_FILE = '2017-03-30_stepFun_Feb28_plateReader_FullFile_1_2017-03-30_11.22.19.mat';
step_down_one = load([STEPDOWN_ONE_DIR '/' STEPDOWN_ONE_FILE]);
step_down_one.up.name = 'step_down_one';
step_down_one.down.name = 'step_down_one';

STEPUPDOWN_TWO_DIR = ['../stepUpDown_repeat_Two_github/saved_data/'];
STEPUPDOWN_TWO_FILE = '2017-03-30_stepFun_repeat2_Mar30_test1_2017-03-30_09.53.38.mat';
step_up_down_two = load([STEPUPDOWN_TWO_DIR '/' STEPUPDOWN_TWO_FILE]);
step_up_down_two.up.name = 'step_up_down_two';
step_up_down_two.down.name = 'step_up_down_two';

%export fig note
exp_fig = @(fH,fName) export_fig([getDate('yyyy-mm-dd') '_' fName '_' ...
    getDate('HH.MM.SS') '.pdf'],...
        '-cmyk','-painters','-pdf', fH);
fncs = {@(y) '', @(y) exp_fig(y{1},y{2})};
expif = @(x,y) feval(fncs{x+1},y);

%% limits between which you want to fit to get l and d coefs (CT hrs)
%loL < l_fit < hiL, loD < l_fit_rng < hiD
%these are in radians. z*CT (CT as in Fig 3CD)
z=((2*pi)/24);
step_up_one.up.lo = z*17; 
step_up_one.up.hi = z*30;
step_up_down_two.up.lo = z*19;
step_up_down_two.up.hi = z*33;

step_down_one.down.lo = z*7;
step_down_one.down.hi = z*21;
step_up_down_two.down.lo = z*7;
step_up_down_two.down.hi= z*21;
clear z;

%% merge expts into a single file 
up_sing = {step_up_one, step_up_down_two};
down_sing = {step_down_one, step_up_down_two};

%adjust phases s.t. ph0=min KaiC phosphorylation, ph_pi = max KaiC
%phosphorylation
for i=1:numel(up_sing)
    phOffset = (7/6)*pi;
    up_sing{i}.up.phase = up_sing{i}.up.phase + phOffset;
    up_sing{i}.up.linPhase = up_sing{i}.up.linPhase + phOffset;
    
    up_sing{i}.up.phase = [up_sing{i}.up.phase-4*pi up_sing{i}.up.phase];
    up_sing{i}.up.phaseShift = [up_sing{i}.up.phaseShift up_sing{i}.up.phaseShift];
end

for i=1:numel(down_sing)
    phOffset = (7/6)*pi;
    down_sing{i}.down.phase = down_sing{i}.down.phase + phOffset;
    down_sing{i}.down.linPhase = down_sing{i}.down.linPhase + phOffset;
    
    down_sing{i}.down.phase = ...
        [down_sing{i}.down.phase-4*pi down_sing{i}.down.phase];
    [~,sortIx] = sort(down_sing{i}.down.phase);
    down_sing{i}.down.phase = down_sing{i}.down.phase(sortIx);
    
    down_sing{i}.down.phaseShift = ...
        [down_sing{i}.down.phaseShift down_sing{i}.down.phaseShift];
    down_sing{i}.down.phaseShift = down_sing{i}.down.phaseShift(sortIx);
end

%mix and match into a 'boot'
for i=1:numel(up_sing)*numel(down_sing)
    pickUp = ceil(i/numel(up_sing));
    pickDown = mod(i,numel(down_sing))+1;
    
    up_mix{i} = up_sing{pickUp}.up;
    down_mix{i} = down_sing{pickDown}.down;
       
    T_hiATP_mix(i) = 0.5*(up_sing{pickUp}.T_hiATP + down_sing{pickDown}.T_hiATP);
    T_loATP_mix(i) = 0.5*(up_sing{pickUp}.T_loATP + down_sing{pickDown}.T_loATP);
end

if TOSAVEMERGE == 1
    save([getDate('yyyy-mm-dd') '_mergedStepFuns_' getDate('HH.MM.SS') '.mat'],...
        'up_mix','down_mix','T_hiATP_mix','T_loATP_mix');
end

%% get wD/wL mean and std
mix_per = {step_up_one, step_down_one, step_up_down_two};
all_hiT = cellfun(@(x) x.T_hiATP, mix_per);
all_loT = cellfun(@(x) x.T_loATP, mix_per);

mu_hi = mean(all_hiT);
mu_lo = mean(all_loT);
sig_hi = std(all_hiT);
sig_lo = std(all_loT);

%wD/wL = T_hi / T_lo
mu_rat = mu_hi / mu_lo;
sig_rat = mu_rat* sqrt((sig_hi/mu_hi)^2 + (sig_lo/mu_lo)^2);

disp(['wD/wL = ' num2str(mu_rat) ' +/- ' num2str(sig_rat)]);

%% get l, d mean and std
all_up = {step_up_one, step_up_down_two};
all_lslope = cellfun(@(x) x.up.linFit(1), all_up); 

all_down = {step_down_one, step_up_down_two};
all_dslope = cellfun(@(x) x.down.linFit(1), all_down);

mu_l = mean(all_lslope);
mu_d = mean(all_dslope);

sig_l = std(all_lslope);
sig_d = std(all_dslope);

disp(['lslope = ' num2str(mu_l) ' +/- ' num2str(sig_l)]);
disp(['dslope = ' num2str(mu_d) ' +/- ' num2str(sig_d)]);

%% re-compute l and d slopes by bootstrapping CT times

%bootstrap l's
for i=1:numel(up_sing)
    stepfun = up_sing{i}.up;
    fitrng = stepfun.phase > stepfun.lo & ...
             stepfun.phase < stepfun.hi;
    x = stepfun.phase(fitrng);
    y = stepfun.phaseShift(fitrng);
    numpts = numel(x);
    nsamples=numpts*100;
    ix = randi(numel(x),nsamples,numpts); 
    for j=1:nsamples      
        sx = x(ix(j,:));
        sy = y(ix(j,:));
        [coef,~] = fitToLine(sx,sy);
        lcoef(j) = coef(1);
    end
end

%bootstrap d's
for i=1:numel(down_sing)
    stepfun = down_sing{i}.down;
    fitrng = stepfun.phase > stepfun.lo & ...
             stepfun.phase < stepfun.hi;
    x = stepfun.phase(fitrng);
    y = stepfun.phaseShift(fitrng);
    numpts = numel(x);
    nsamples=numpts*100;
    ix = randi(numel(x),nsamples,numpts); 
    for j=1:nsamples      
        sx = x(ix(j,:));
        sy = y(ix(j,:));
        [coef,~] = fitToLine(sx,sy);
        dcoef(j) = coef(1);
    end
end

mslopefun = (@(L,D,OMR) (1- (-(1-L)*OMR/(D+L-L*D)) - (1/(D+L-L*D))));

%bootstrap mslope
OMR=1; %omegaD/omegaL
for j=1:nsamples
    mslope(j) = mslopefun(lcoef(j),dcoef(j),OMR);
end

mu_l = mean(lcoef);
sig_l = std(lcoef);
mu_d = mean(dcoef);
sig_d = std(dcoef);
mu_m = mean(mslope);
sig_m = std(mslope);
disp(['boot lslope = ' num2str(mu_l) ' +/- ' num2str(sig_l)]);
disp(['boot dslope = ' num2str(mu_d) ' +/- ' num2str(sig_d)]);
disp(['boot mslope = ' num2str(mu_m) ' +/- ' num2str(sig_m) ...
    ', for OMR=' num2str(OMR)]);


%break;

%% make a double suplot figure comparing the old L and D and the new D

z=24/(2*pi);

%
fCompUp = figure();
colors = {'b','c'};
dsUp=[];
for i=1:numel(up_sing)
    ptstosave = (z*up_sing{i}.up.phase) > -1;
    numpts = 1:sum(ptstosave);
    dsUp(numpts,2*(i-1)+1) = ...
        (z*up_sing{i}.up.phase(ptstosave))';
    dsUp(numpts,2*(i-1)+2) = ...
       (z*up_sing{i}.up.phaseShift(ptstosave))';
   
    p(i) = plot(z*up_sing{i}.up.phase, z*up_sing{i}.up.phaseShift, 's-',...
        'color',colors{i},'markersize',4,'markerfacecolor',colors{i});
    hold on;
%     pLin(i) = plot(z*up_sing{i}.up.linPhase,z*up_sing{i}.up.linPhaseShift','-',...
%         'color',colors{i},'marker','none');
end
dsUp = mat2dataset(dsUp);
dsUp.Properties.VarNames = {'Rep1_stepTime','Rep1_phaseShift',...
    'Rep2_stepTime','Rep2_phaseShift'};
export(dsUp,'File','stepUP_polarization.csv','Delimiter',',');

set(gca,'xtick',z*[-6*pi:0.5*pi:6*pi],'xlim',z*([-2*pi 2*pi]+2*pi),...
    'ylim',z*[-pi pi], 'ytick', z*[-2*pi:(1/3)*pi:2*pi]);
grid off;
xlabel('step time (CT hours)');
ylabel('phase shift (CT hours)');

fCompDown = figure();
colors = {'r','m'};
dsDown = [];
for i=1:numel(down_sing)
    ptstosave = (z*down_sing{i}.down.phase) > -1;
    numpts = 1:sum(ptstosave);
    dsDown(numpts,2*(i-1)+1) = ...
        (z*down_sing{i}.down.phase(ptstosave))';
    dsDown(numpts,2*(i-1)+2) = ...
       (z*down_sing{i}.down.phaseShift(ptstosave))';
    p(i) = plot(z*down_sing{i}.down.phase, z*down_sing{i}.down.phaseShift, 's-',...
        'color',colors{i},'markersize',4,'markerfacecolor',colors{i});
    hold on;
%     pLin(i) = plot(z*down_sing{i}.down.linPhase,z*down_sing{i}.down.linPhaseShift','-',...
%         'color',colors{i},'marker','none');
end
dsDown = mat2dataset(dsDown);
dsDown.Properties.VarNames = {'Rep1_stepTime','Rep1_phaseShift',...
    'Rep2_stepTime','Rep2_phaseShift'};
export(dsDown,'File','stepDown_polarization.csv','Delimiter',',');

set(gca,'xtick',z*[-6*pi:0.5*pi:6*pi],'xlim',z*([-2*pi 2*pi]+2*pi),...
    'ylim',z*[-pi pi], 'ytick', z*[-2*pi:(1/3)*pi:2*pi]);
grid off;
xlabel('step time (CT hours)');
ylabel('phase shift (CT hours)');

set(fCompUp,'units','inches','position',[3 3 3.75 3]);
set(fCompDown,'units','inches','position',[3 3 3.75 3]);


%% load metaboloperiod data
JOINT_FIT_DIR = ['.'];
JOINT_FIT_FILE = ['2017-03-17_inVitroPP_fits_eachExpOwnPer_03.53.51.mat'];
vitro = load([JOINT_FIT_DIR '/' JOINT_FIT_FILE], 'vitro');
vitroPR=vitro.vitro;

%%
good = [1:5 12]; %pick pp=6-12 from first expt, pp=16 from second expt
pp = vitroPR.pp(good);
colors = hot(numel(pp));

%recall that both step funs and fit data is in radians
dawn = vitroPR.dawnPhase + phOffset; %this is really theta(dawn)+L(theta(dawn))
dawn = mod(dawn(good),2*pi);
dawn = wrapVecAround(dawn,1.5*pi,2*pi,'gt');
dawn = wrapVecAround(dawn,pi,2*pi,'lt')

%periods in light and dark
TLIGHT = up_sing{1}.T_hiATP;
TDARK = up_sing{1}.T_loATP;

pickD = 1; %which D fun to use?

%compute times of dawn and dusk to be marked with arrows
%dusk = dawn + pp*1/TLIGHT
dusk = dawn + (2*pi/TLIGHT)*pp;
dusk = wrapVecAround(dusk,0,2*pi,'lt');
dusk = wrapVecAround(dusk,8,2*pi,'gt');
dawnStep = interp1(up_sing{1}.up.phase, up_sing{1}.up.phaseShift, dawn);
duskStep = interp1(down_sing{pickD}.down.phase,...
    down_sing{pickD}.down.phaseShift,dusk);

%mark colors
for i=1:numel(dusk)
    figure(fCompDown);
    z=24/(2*pi); %this z must be in the same units as in the plots above
    
    %try arrows
    arrow(z*[dusk(i) -0.5*pi],z*[dusk(i) -0.25*pi],'length',8,...
        'facecolor',colors(i,:),'edgecolor','k','tipangle',20);
    hold on;
    
    figure(fCompUp);
    %try arrows
    arrow(z*[dawn(i) -0.5*pi],z*[dawn(i) -0.25*pi],'length',8,...
        'facecolor',colors(i,:),'edgecolor','k','tipangle',20);
    hold on;
    
end

%save arrows as datasets
dsDusk = [pp' z*dusk' z*-0.25*pi*ones(size(dusk'))];
dsDusk_names = {'day_length','dusk_time_CT', 'dusk_arrow_top_CT'};
dsDusk = mat2dataset(dsDusk);
dsDusk.Properties.VarNames = dsDusk_names;
export(dsDusk,'FILE','dusk_arrows.csv','Delimiter',',');

dsDawn = [pp' z*dawn' z*-0.25*pi*ones(size(dusk'))];
dsDawn_names = {'day_length','dawn_time_CT', 'dawn_arrow_top_CT'};
dsDawn = mat2dataset(dsDawn);
dsDawn.Properties.VarNames = dsDawn_names;
export(dsDawn,'FILE','dawn_arrows.csv','Delimiter',',');

%export figures?
expif(TOEXP,{fCompUp,'Lfun'});
expif(TOEXP,{fCompDown,'Dfun'});
