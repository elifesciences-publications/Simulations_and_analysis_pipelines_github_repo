%2016-08-31, EL: modified to accept L and D funs as inputs. Pass them to
%getInVitroPhaseShiftBoot.m directly as STEPUP and STEPDOWN. STEPUP is a
%struct with fields STEPUP.phase, STEPUP.phaseShift. Likewise for STEPDOWN.
%2016-08-28, EL: modified to allow for bootstrapped L and D funs
%2016-08-01, EL: modified to make plots for manuscript
%Drive a phase oscillator with 'light/dark' cycles according to
%stepup/stepdown functions (see code re: response to dawn and dusk). Input
%period of oscillator in light/dark, driving period, fraction of time
%lights are on and the duration of simulation.
%   Input: TDRIVE, TLIGHT, TDARK in hours. DUTYFRAC in [0,1], TEND in hrs.
%   Output: DAWNPHASE [0,1] of last dawn before TEND and corresponding
%   DAWNPHASESHIFT in response to the last dawn. Likewise for DUSK.
function [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
    drivePhaseOscilStepFn_WithBoot(TDRIVE, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
    DUTYFRAC, TEND, NUMCYC, STARTPHASE)
%% TEST
TOTEST=0;
if TOTEST==1
    fname = mfilename('fullpath');
    warning([fname ': Running in test mode!']);
    
    STEPFUNTYPE='BOOT';
    switch STEPFUNTYPE
        case 'KaiTEST'
            %KaiTEST values
            TDRIVE=26.7;
            TLIGHT=26.7;
            TDARK=20.94;
        case 'InVitro'
            %in vitro values
            TLIGHT=23.60;
            TDARK=22.51;
            TDRIVE=24.0;
        case 'BOOT'
            INDIR='/Users/E/Documents/Advisers/Rust/Clocks/2015-11-23_stepUp-stepDown/2016-10-31_stepFunOct';
            INFILE='2016-10-31_stepFun_Lin-NLin-Oct_1000_2016-10-31_20.45.44.mat'; %700 boots
            load([INDIR '/' INFILE]);
            BOOTNUM=111;
            TLIGHT=T_hiATP(BOOTNUM);
            TDARK=T_loATP(BOOTNUM);
            
            STEPUP.phase = up.phase(BOOTNUM,:);
            STEPUP.phaseShift = up.phaseShift(BOOTNUM,:);
            STEPDOWN.phase = down.phase(BOOTNUM,:);
            STEPDOWN.phaseShift = down.phaseShift(BOOTNUM,:);
            
            TDRIVE=24.0; 
        otherwise
            error([fname ': Did not decide on a valid StepFun type!']);
    end
    
    %driving params to both
    DUTYFRAC=[12]/24;
    NUMCYC=1:20;
    TEND=TDRIVE*max(NUMCYC)+1;
    STARTPHASE=0;
end

%% make phase oscillator predictions using interpolated stepup/down functions
TODISP=1; %suppress print output during simulations

Tlight = TLIGHT; %upendper (hrs);
Tdark = TDARK; %upstartper (hrs);

Tdrive=TDRIVE;
dutyDrive=DUTYFRAC;
simtime = TEND; %hrs

%get step func. handle
KAITEST_STEPFUN='getPhaseShift';
INVITRO_STEPFUN='getInVitroPhaseShift';
BOOT_STEPFUN = 'getInVitroPhaseShiftBoot';
stepFun = str2func(BOOT_STEPFUN);
dispif(TODISP,['Using ' func2str(stepFun) ' in ' mfilename('fullpath') '.']);

%start frpm ph[0]=0 by default
if isempty(STARTPHASE)
    ph=[0];
else
    ph=[STARTPHASE];
    dispif(TODISP,['ph0=' num2str(STARTPHASE)]);
    if STARTPHASE > 1 | STARTPHASE < 0
        warning(['Starting phase = ' num2str(STARTPHASE) ' is outside [0,1]!']);
    end
end

dayInd=[1];
t=[];
lightson=[];
dusk=[];
dawn=[];

dt = 0.01; %must set to <<0.1 because when driving with fractions of 26.7
            %need to be wary of rounding error! May skip dawn or dusk!

t=0:dt:simtime;
for tind=1:numel(t)
    islight = (mod(t(tind),Tdrive) <= dutyDrive*Tdrive+0.5*dt &...
        mod(t(tind),Tdrive) >= 0.5*dt); %dusk is last light pt; dawn is last dark pt
    isdawn = mod(t(tind),Tdrive) >= 0-0.5*dt & mod(t(tind),Tdrive) < 0.5*dt;
    isdusk = mod(t(tind),Tdrive) > dutyDrive*Tdrive-0.5*dt & ...
        mod(t(tind),Tdrive) <= dutyDrive*Tdrive+0.5*dt;
%     disp(['t=' num2str(t(tind)) ...
%         '; isdusk=' num2str(isdusk) ...
%         '; isdawn=' num2str(isdawn)]);
    
    if islight || tind==1
        ph(tind+1) = ph(tind) + (1/Tlight)*dt;
        lightson(tind) = 1;
    else
        ph(tind+1) = ph(tind) + (1/Tdark)*dt;
        lightson(tind) = 0;
    end
    
    dawn(tind+1) = 0;
    dawnshift(tind+1) = 0;
    dusk(tind+1) = 0;
    duskshift(tind+1) = 0;
    
    if tind < numel(t)
        dayInd(tind+1) = dayInd(tind);
    end

    if isdawn & tind>1 & lightson(tind-1) == 0 % must transition from dark to light
        %disp(['dawn at t=' num2str(t(tind)) ',ph=' num2str(ph(tind+1))]);
        %disp(['dawn pshift=' num2str(getPhaseShift(mod(ph(tind+1),1),'stepup'))]);
        dawn(tind) = 1;
        %dawnshift(tind) = getPhaseShift(mod(ph(tind),1),'stepup');
        dawnshift(tind) = stepFun(mod(ph(tind),1),STEPUP);
        ph(tind+1) = ph(tind+1) + dawnshift(tind);
        dayInd(tind+1) = dayInd(tind)+1;
       
    elseif isdusk & lightson(tind-1) == 1
        %disp(['dusk at t=' num2str(t(tind))]);
        dusk(tind) = 1;
        %duskshift(tind) = getPhaseShift(mod(ph(tind),1),'stepdown');
        duskshift(tind) = stepFun(mod(ph(tind),1),STEPDOWN);
        ph(tind+1) = ph(tind+1) + duskshift(tind);
    end
    
end

t(end+1) = t(end)+dt;
lightson=logical(lightson);

%sanity check that haven't skipped a dawn or dusk
if ((sum(dusk) ~= sum(dawn)) || (sum(dusk) ~= max(NUMCYC))) && (TODISP == 1)
    fname = mfilename('fullpath');
    warning([fname ': PO missed a dusk or dawn!']);
end

%find all pts right after lights on, include stepup response
% daystart=find(logical(dawn))+1; 
% DAWNPHASE = ph(daystart);
if DUTYFRAC*Tdrive > 0+dt & DUTYFRAC*Tdrive < Tdrive-dt
    DAWNPHASE = ph(logical(dawn));
    DAWNPHASE = mod(DAWNPHASE(NUMCYC),1);
    DAWNPHASESHIFT = dawnshift(logical(dawn)); %pshift at last dawn/dusk?
    DAWNPHASESHIFT = DAWNPHASESHIFT(NUMCYC);
end

% nightstart=find(logical(dusk))+1;
% DUSKPHASE = ph(nightstart);
if DUTYFRAC*Tdrive > 0+dt & DUTYFRAC*Tdrive < Tdrive-dt
    DUSKPHASE = ph(logical(dusk));
    DUSKPHASE = mod(DUSKPHASE(NUMCYC),1);
    DUSKPHASESHIFT = duskshift(logical(dusk));
    DUSKPHASESHIFT = DUSKPHASESHIFT(NUMCYC);
end

TO_EXPORT_LDCYCLES_FIG = 0;
TOMAKELDCYCLES_FIG = 0;
%plot light-dark boxes overlaid on top of phase oscillator trajectories
if TOMAKELDCYCLES_FIG == 1
    fLDTraj = figure();
    
    for n=NUMCYC
        todayLight = lightson & (dayInd == n);
        todayDark = (~lightson) & (dayInd == n);
        pLight(n) = plot(t(todayLight),mod(ph(todayLight),1),'g.','markersize',6,'linewidth',1);
        hold on;
        pDark(n) = plot(t(todayDark),mod(ph(todayDark),1),'g.','markersize',6,'linewidth',1);

%         %if want to plot this in sine coordinates instead 
%         pLight(n) = plot(t(todayLight),sin(2*pi*(ph(todayLight)-0.5)),'r','linewidth',2);
%         hold on;
%         pDark(n) = plot(t(todayDark),sin(2*pi*(ph(todayDark)-0.5)),'r','linewidth',2);

        greyCol = [209 210 212]/255;
        yellowCol = [235 235 235]/255;
        
        pDarkBox = area(t(todayDark),...
            ones(size(t(todayDark))),...
            'facecolor',greyCol,'edgecolor','none');
               
        if n ~=1
            pLightBox = area(t(todayLight),...
                ones(size(t(todayLight))),...
                'facecolor',yellowCol,'edgecolor','none');
        else
            
            pLightBox = area([0 t(todayLight)],...
                [1 ones(size(t(todayLight)))],...
                'facecolor',yellowCol,'edgecolor','none');
        end
        uistack(pDark(n),'top');
        uistack(pLight(n),'top');
    end
    
    %picture size [5x3] for illustrator 3 panels
    set(gca,'color','none','layer','top'); %to make sure axes are above colored area bars
    set(gca,'xtick',0:12:240,'ytick',[0:0.25:1],'xticklabel',[0:12:240],'fontsize',12);
    xlabel('time (hours)');
    ylabel('oscillator phase \theta (rad/2\pi)');
    grid off;
    %xlim([0 72]);
    set(gcf,'units','inches','position',[0 0 8.5 2.5],'color','none');
    
    if TO_EXPORT_LDCYCLES_FIG == 1
        savename=[getDate('yyyy-mm-dd') '_PhOsc_Illustration_numcyc-' ...
            num2str(max(NUMCYC)) '_df-' num2str(DUTYFRAC) '_' ...
            datestr(now, 'yyyy-mm-dd_HH.MM.SS') ...
            '.pdf'];
        export_fig(savename,'-cmyk','-painters',fLDTraj);
    end
end

%plot phase oscillator trajectories, mark dawn and dusk
TOPLOT=0;
if TOPLOT==1
    fPhaseOsc = figure();
    plot(t(lightson),mod(ph(lightson),1),'r.','markersize',7,'linewidth',1);
    hold on;
    plot(t(~lightson),mod(ph(~lightson),1),'k.','linewidth',1);
    hold on;
    
    if DUTYFRAC*Tdrive > 0+dt & DUTYFRAC*Tdrive < Tdrive-dt
        %dusk
        plot(t(logical(dusk)),mod(ph(logical(dusk)),1),...
            'mo','markersize',7,'linewidth',1);
        hold on;
        
        %dawn
        plot(t(logical(dawn)),mod(ph(logical(dawn)),1),...
            'go','markersize',7','linewidth',1);
        hold on;
        
        %dawn+dawnshift
        plot(t(logical(dawn)),...
            mod(...
            ph(logical(dawn))+dawnshift(logical(dawn))+(1/Tlight)*dt,...
            1),...
            'gs','markerfacecolor','g','markersize',7','linewidth',1);
        hold on;
        
        %dusk+duskshift
        plot(t(logical(dusk)),...
            mod(...
            ph(logical(dusk))+duskshift(logical(dusk))+(1/Tdark)*dt,...
            1),...
            'ms','markerfacecolor','m','markersize',7','linewidth',1);
        hold on;
    end

    titStr = ['Tlight=' num2str(Tlight) ', Tdark=' num2str(Tdark) ...
        10 'Tdrive=' num2str(Tdrive) ', dutyDrive=' num2str(dutyDrive)]; 
    
    if DUTYFRAC*Tdrive > 0+dt & DUTYFRAC*Tdrive < Tdrive-dt
        titStr=[titStr ...
        10 'end dawnPhase=' sprintf('%2.4f',DAWNPHASE(end)) ',' ...
        10 'end dawnPshift=' sprintf('%2.4f',DAWNPHASESHIFT(end))];
    end
    
    title(titStr);
    xlim([min(t) max(t)]);
    set(gca,'xtick',min(t):24:max(t),'ytick',[0:0.25:1]);
    xlabel('hours');
    ylabel('clock phase (rad/2pi)');    
  
end