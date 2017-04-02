%2017-03-28, EL: global fit of datasets in Fig. 5
%dependencies: export_fig.m from FileExchange to export figures

clear all;
close all;

cd('.');

%% load files 
load('2016-10-25_16.28.08_prcdata.mat'); %[dpDur dpTime pshift pshiftErr pshiftStdev ct_time norm_phase]
load('2016-10-25_16.56.02_TeeTauData.mat'); %[tau T pkLL stderrPkLL stdevPkLL]
load('2016-10-25_23.24.46_wedgedata.mat');
ppdata = TeeTauData;
clear TeeTauData;

%% switch order of columns in wedge and prc
prcdata = prcdata(:, [2 1 3:end]);
wdata = wdata(:,[2 1 3:end]);
wdata(:,1) = wdata(:,1) - 24;
%% set up fit funs

%entrained phase: 
%X(:,1) = tau, X(:,2) = T; X(:,3) = peakTime in LL. 
%peakTime(tau,T) = tau*(1-b1-b2) + T*b1 + b3;
peakTimeFun = (@(b,X) X(:,1)*(1-b(1)-b(2)) + X(:,2)*b(1) + b(3));

%phase shift:
%X(:,1) = CT, X(:,2) = delta;
%pshift(CT,delta) = CT*(1/b2) + delta*(1+b1/b2) + b4; (CT = theta_t/omega_L)
phaseShiftFun  = (@(b,X) X(:,1)*(1/b(2)) + X(:,2)*(1+(b(1)/b(2))) + b(4));
phaseShiftFun2 = (@(b,X) X(:,1)*(1/b(2)) + X(:,2)*(1+(b(1)/b(2))) + b(5));
phaseShiftFun3 = (@(b,X) X(:,1)*(1/b(2)) + X(:,2)*(1+(b(1)/b(2))) + b(6));

%wedge alone:
%X(:,1) = CT; X(:,2) = delta;
%wedgeshift(delta) = delta*(1+b1/b2) + b4
% wedgeFun = (@(b,X) X(:,2)*(1+(b(1)/b(2))) + b(4));
% wedgeFun2 = (@(b,X) X(:,2)*(1+(b(1)/b(2))) + b(5));
wedgeFun  = phaseShiftFun3;
wedgeFun2 = phaseShiftFun3;

%prc alone:
%X(:,1) = CT;
%prcshift(ct) = ct*(1/b2) + b5;
% prcFun = (@(b,X) X(:,1)*(1/b(2)) + b(6));
% prcFun2 = (@(b,X) X(:,1)*(1/b(2)) + b(7));
prcFun = phaseShiftFun;
prcFun2 = phaseShiftFun2;

%% call fitter

%set up data structures for fitter
data.pp = ppdata;
fitfun.pp = peakTimeFun;
%
data.wedge = wdata(wdata(:,1) == 60-24,:);
fitfun.wedge = wedgeFun;
% 
data.wedge2 = wdata(wdata(:,1) == 63-24,:);
fitfun.wedge2 = wedgeFun2;
%
prcdata(prcdata(:,2) < -6,2) = prcdata(prcdata(:,2) < -6,2) + 24;
prcdata(prcdata(:,2) > 12,2) = prcdata(prcdata(:,2) > 12,2) - 24;
prcdata(:,1) = prcdata(:,1) - 24;
pickprc=prcdata(:,1) > 33 & prcdata(:,1) < 46;
pickprc2=prcdata(:,1) > 51 & prcdata(:,1) < 62;

data.prc = prcdata(pickprc,:);
fitfun.prc = prcFun;
data.prc2 = prcdata(pickprc2,:);
fitfun.prc2 = prcFun2;
%

%set up guesses for regression params
%beta0 = [-1.4288    1.9269   42.6946] optimized for beta1-beta3 for T and tau
%beta0 = [-1.4211    1.9261   42.5857   -4.7074   -2.4109] optimized for beta1-beta5 for T,
%tau and wedge
%beta0=[-1.2823    1.7449   38.7851   -4.7432   -2.4467   -6.7788] for
%T,tau, wedge and PRC -- separate y intercepts
beta0 = [-1.2823    1.7449   38.7851   -4.7432   -2.4467   -6.7788];
lb = -Inf*ones(size(beta0));
ub = Inf*ones(size(beta0));
options = optimoptions('lsqnonlin','ScaleProblem','Jacobian');
options.MaxFunEvals = 1000;

%% fit
[beta_fit, resnorm, residual, exitflag, output ,lambda, jacobian] = ...
    lsqnonlin(@globalMultiLinFitFun,beta0,lb,ub,...
    options,...
    data, fitfun);

%% get errors on parameter estimates
dof = numel(residual) - numel(beta_fit);
chisq = resnorm/dof;
chisqProb = chi2cdf(resnorm,dof); %prob. chisqRV < chisq

j=full(jacobian);
varpar =chisq*(inv(j'*j));
sigma = sqrt(diag(varpar));

confint = nlparci(beta_fit,residual,'jacobian',jacobian);

%get variables for saving
out.BETA_FIT = beta_fit;
out.BETA_SIGMA = sigma;
out.CHISQ = chisq;
out.CHISQ_PROB = chisqProb;
out.CI = confint; %checked that confint agrees with sigma

disp(['tauSlope = ' num2str((1-beta_fit(1)-beta_fit(2))) ...
    ' +/- ' num2str(sqrt(sigma(1)^2 + sigma(2)^2))]);
disp(['TSlope = ' num2str(beta_fit(1)) ...
    ' +/- ' num2str(sqrt(sigma(1)^2))]);

%% test fit

%use stderr or stdev in plots
sigma=5;  %5 = stdev, 4 = stderr

%for phase vs pp
fPeakTime = figure();
uniT = unique(ppdata(:,2));
colors = hsv(numel(uniT));
pickUniT=[1:5];
for u=pickUniT
    subset = ppdata(:,2) == uniT(u);
    plterr(u) = errorbar(ppdata(subset,1), ppdata(subset,3), ...
        ppdata(subset,sigma),'s','color',colors(u,:),...
        'markerfacecolor',colors(u,:),...
        'markersize',4,'markeredgecolor','none','linewidth',1);
    hold on;

    pltfit(u) = plot(ppdata(subset,1),peakTimeFun(beta_fit,ppdata(subset,:)),...
        'marker','none','color',colors(u,:),'linewidth',2);
    lab{u} = ['T=' num2str(uniT(u)) ' hours'];
end
xlabel('day length \tau (hours)');
ylabel('peak time t_{pk} (hours)');
set(gca,'xlim',[6 18],'ylim',[6 18],'fontsize',12);
legend(plterr(pickUniT),lab{pickUniT},'location','southeast');
legend boxoff;
grid off;

%for wedge exp
fWedge = figure();
uniPh = unique(wdata(:,1));

mySet = wdata(:,1) == uniPh(1);
wplt(1) = errorbar(wdata(mySet,2), wdata(mySet,3), wdata(mySet,sigma),...
    's','markerfacecolor','b','color','b',...
    'markersize',4,'markeredgecolor','none','linewidth',1);
hold on;
wfitplt1 = plot(wdata(mySet,2),wedgeFun(beta_fit,wdata(mySet,:)),'b','linewidth',2);

mySet = wdata(:,1) == uniPh(2);
wplt(2) = errorbar(wdata(mySet,2), wdata(mySet,3), wdata(mySet,sigma),...
    '>','markerfacecolor','b','color','b',...
    'markersize',4,'markeredgecolor','none','linewidth',1);
wfitplt2 = plot(wdata(mySet,2),wedgeFun(beta_fit,wdata(mySet,:)),'b','linewidth',2);

xlabel(['dark pulse length '  '\delta (hours)']);
ylabel(['phase shift '  '\Deltat_{pk}(hours)']);
set(gca,'xlim',[6 18],'fontsize',12);
legend(wplt,{...
    ['t_{DP}=' num2str(wdata(1,1)) ' hours ' ...
    '(\theta=' num2str(wdata(1,end),'%2.2f') ')'], ...
            ['t_{DP}=' num2str(wdata(end,1)) ' hours (\theta=' ...
            num2str(wdata(end,end),'%2.2f') ')']},...
            'location','northwest');
legend boxoff;
box on;
grid off;
     
%for prc expt
fPRC = figure();
%
prcplterr = errorbar(prcdata(:,1),prcdata(:,3),prcdata(:,sigma),'s',...
    'markerfacecolor','k','color','k',...
    'markersize',4,'markeredgecolor','none','linewidth',1);
hold on;
%
mySet = pickprc;
prcfitplt(1) = plot(prcdata(mySet,1),prcFun(beta_fit,prcdata(mySet,:)),'-',...
    'linewidth',2,'color','k');
mySet = pickprc2;
prcfitplt(2) = plot(prcdata(mySet,1),prcFun2(beta_fit,prcdata(mySet,:)),'-',...
    'linewidth',2,'color','k');
%
set(gca,'xtick',33:3:63,'ytick',-12:4:15,'fontsize',12);
legPRC=legend(['\delta=12 hours'],'location','northwest');
legend boxoff;
xlim([36-2 60+2]);
ylim([-10-2 14+2]);
xlabel(['dark pulse time '  't_{DP} (hours)']);
ylabel(['phase shift ' '\Deltat_{pk}(hours)']);
grid off;

%% save
TOSAVE = 0;
if TOSAVE == 1
formatFig = (@(f) set(f,'units','inches','position',[0 0 3.75 3.75]));
expFig = (@(f,fname) export_fig(fname, '-cmyk','-painters','-pdf',f));
getDate = @(x) datestr(now,'yyyy-mm-dd_HH.MM.SS');

formatFig(fPeakTime);
formatFig(fWedge);
formatFig(fPRC);

expFig(fPeakTime,['tauTFig-all_stderr' getDate()]);
expFig(fWedge,['wedgeFig-allBlue_stderr' getDate()]);
expFig(fPRC,['prcFig-allBlack_stderr' getDate()]);
save(['fitParams_' getDate() '.mat'], 'out');

prclab = [prcdata(:,1) prcdata(:,end)];
%save(['prcPhaseLabels_' getDate() '.mat'], 'prclab');

end
