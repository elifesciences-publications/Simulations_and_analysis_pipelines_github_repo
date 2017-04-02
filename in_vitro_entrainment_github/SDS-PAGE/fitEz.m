%2016-03-30, EL: produces fits for manuscript Fig 2C.

clear all;
close all;
clc;

%% load data 
DATA = loadQuants('in_vitro_entrainment_SDS-PAGE_raw_data.csv');
time = DATA(:,3);
DATASHORT = DATA(time > 48,:); %fit only day 3
rxnset = DATASHORT(:,1);

%% fit to sine function with fixed 24-hr period
SINEFUN = str2func('ezSine_NoPer');

%% fit
%no period [offset ampl phase]
beta0 = [50 50 pi]';
lb = [0 0 0]';
ub = [100 100 2*pi]';
phaseIx = 3;

beta=[];

uniset = unique(rxnset);
for r=1:numel(uniset)
    thisset = DATASHORT(:,1) == uniset(r);
    x{r}=DATASHORT(thisset,3);
    y{r}=DATASHORT(thisset,end);
       
    %do first fit to generate guesses for second iteration
    [guessfit, ~, ~, ~, ~, ~, ~] = lsqnonlin(SINEFUN,beta0,lb,ub,[],{x{r}},{y{r}});
    
    guesslb = lb;
    guesslb(phaseIx) = guessfit(phaseIx) - pi;
    guessub = ub;
    guessub(phaseIx) = guessfit(phaseIx) + pi;
    
    %fit second time using updated initial guess and phase limits
    options = optimoptions('lsqnonlin','Jacobian','on');
    %options = [];

    [fit{r},resnorm{r},residual{r},exitflag{r},output{r},lambda{r},jacobian{r}] =...
        lsqnonlin(SINEFUN,guessfit,guesslb,guessub,options,{x{r}},{y{r}});
    fival{r} = SINEFUN(fit{r},{x{r}},0);

    ph(r) = fit{r}(phaseIx); %mod((0.5*pi+fit{r}(4))*(12/pi),24);
    
    %this photoperiod
    thesepp = DATASHORT(thisset,2);
    pp_plot(r) = thesepp(1);
    
    %use nlparci for confidence intervals
    ci{r} = nlparci(fit{r},residual{r},'jacobian',jacobian{r},'alpha',0.33);
    lo_ph(r) = ci{r}(phaseIx,1); %(0.5*pi+ci{r}(4,1))*(12/pi);
    hi_ph(r) = ci{r}(phaseIx,2); %(0.5*pi+ci{r}(4,2))*(12/pi);
    
    %compute stdev of phase by hand
    dof{r} = numel(residual{r}) - numel(fit{r});
    varp{r} = (resnorm{r}/dof{r})*(inv(jacobian{r}'*jacobian{r}));
    sigma{r} = sqrt(diag(varp{r}));
    ph_err(r) = sigma{r}(phaseIx); %(12/pi)*sigma{r}(4);
end

for r=1:numel(uniset)
    subplot(3,3,r);
    plot(x{r},y{r},'bo');
    hold on;
    plot(x{r},fival{r},'b-');
end

%% plot fits

%do modular arithmetic
ph = mod(ph, 2*pi);
ph = wrapVecAround(ph,0.75*2*pi, 2*pi, 'gt');

lo_ph = mod(lo_ph,2*pi);
hi_ph = mod(hi_ph, 2*pi);

for n=1:numel(ph)
    lo_ph(n) = wrapVecAround(lo_ph(n),ph(n),2*pi,'gt');
    hi_ph(n) = wrapVecAround(hi_ph(n), ph(n), 2*pi, 'lt');
end

%fit to line
fLinFit = figure();
z=24/(2*pi);
goodpp = pp_plot < 24;
ph_err = full(ph_err);
[lincoef, errcoef] = ...
    fitToLineYErr(pp_plot(goodpp),z*ph(goodpp), z*ph_err(goodpp));

errorbar(pp_plot(goodpp),z*ph(goodpp), z*ph_err(goodpp), 'k');
hold on;
plot(pp_plot(goodpp),polyval(lincoef,pp_plot(goodpp),'b-'));

%% collect and save data
vitro.pp = pp_plot;
vitro.fitPhase = ph;
vitro.fitPhase_Err = ph_err;
vitro.dawnPhase_relMin = (-ph-0.5*pi);
vitro.dawnPhase_01 = vitro.dawnPhase_relMin/(2*pi);
vitro.dawnPhase_01_Err = vitro.fitPhase_Err/(2*pi);
vitro.pkT = (24/(2*pi))*(0.5*pi-(-vitro.fitPhase));
vitro.pkT_Err = (24/(2*pi))*vitro.fitPhase_Err;
vitro.good = goodpp;

%do fits
z=(24/(2*pi));
[lincoef, errcoef] = ...
    fitToLineYErr(vitro.pp(goodpp),...
                    vitro.pkT(goodpp), vitro.pkT_Err(goodpp));
vitro.pkT_lincoef = lincoef;
vitro.pkT_errcoef = errcoef;

[lincoef_ph, errcoef_ph] = ...
    fitToLineYErr(vitro.pp(goodpp),...
                    vitro.dawnPhase_01(goodpp), vitro.dawnPhase_01_Err(goodpp));
vitro.ph_lincoef = lincoef_ph;
vitro.ph_errcoef = errcoef_ph;

%plot
figure();
subplot(2,1,1);
errorbar(vitro.pp,vitro.dawnPhase_01,vitro.dawnPhase_01_Err,'o');
hold on;
plot(vitro.pp(goodpp),polyval(vitro.ph_lincoef,vitro.pp(goodpp)));
title('dawn phase');
subplot(2,1,2);
errorbar(vitro.pp, vitro.pkT, vitro.pkT_Err,'o');
hold on;
plot(vitro.pp(goodpp), polyval(vitro.pkT_lincoef,vitro.pp(goodpp)));
title('pk time');

%save
%save([getDate('yyyy-mm-dd') '_vitro_day3_' getDate() '.mat'],'vitro');