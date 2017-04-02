%2017-03-28, EL: run phase-resetting simulations for a phase oscillator
%You can run this code either as a function (called by
%drivePhaseOscilStepFn_PRC.m or as a script, by setting TOTEST=1 and
%manually inputting parameters you want to check in your simulation or to
%generate an illustration (such as in Math Appendix).

%Drive a phase oscillator with 'light/dark' cycles according to
%stepup/stepdown functions (see code re: response to dawn and dusk). Input
%period of oscillator in light/dark, driving period, fraction of time
%lights are on and the duration of simulation.
%   Input: TDRIVE, TLIGHT, TDARK in hours. DUTYFRAC in [0,1], TEND in hrs.
%   Output: DAWNPHASE [0,1] of last dawn before TEND and corresponding
%   DAWNPHASESHIFT in response to the last dawn. Likewise for DUSK.
function [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
    drivePhaseOscilStepFn_PRC(TEND, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
    DPTIME, DPDUR, STARTPHASE)
%% TEST
TOTEST=0;
if TOTEST==1
    fname = mfilename('fullpath');
    warning([fname ': Running in test mode!']);
    
    %input file to load step functions       
    INDIR=['../helper functions and shared files/'];
    INFILE = '2017-04-01_mergedStepFuns_11.35.34.mat.mat';
    load([INDIR '/' INFILE]);
    
    %if your step function file contains multiple step functions (e.g.,
    %bootstraps (%P KaiC datasets) or multiple L and D combinations
    %(polarization datasets)), select which one you want to use here
    BOOTNUM=1; %all examples in manuscript ran using BOOTNUM=1
    TLIGHT=T_hiATP_mix(BOOTNUM);
    TDARK=T_loATP_mix(BOOTNUM);
    STEPUP.phase = up_mix{BOOTNUM}.phase;
    STEPUP.phaseShift = up_mix{BOOTNUM}.phaseShift;
    STEPDOWN.phase = down_mix{BOOTNUM}.phase;
    STEPDOWN.phaseShift = down_mix{BOOTNUM}.phaseShift;
    
    TDRIVE=24.0;
    
    %driving params to both
    DPTIME=7; %hrs
    DPDUR=9; %hrs
    TEND=TDRIVE*3;
    STARTPHASE=0;
end

%% export settings
TO_EXPORT_LDCYCLES_FIG = 0;
TOMAKELDCYCLES_FIG = 0; %do you want to make the figure illustrating phase resetting?

%% make phase oscillator predictions using interpolated stepup/down functions
Tlight = TLIGHT; %upendper (hrs);
Tdark = TDARK; %upstartper (hrs);

dptime = DPTIME;
dpdur = DPDUR;
simtime = TEND; %hrs

%get step func. handle
BOOT_STEPFUN = 'getInVitroPhaseShiftBoot';
stepFun = str2func(BOOT_STEPFUN);
disp(['Using ' func2str(stepFun) ' in ' mfilename('fullpath') '.']);

%start frpm ph[0]=0 by default
if isempty(STARTPHASE)
    ph=[0];
    phControl=[0]; %for LL curve
else
    ph=[STARTPHASE];
    phControl=[STARTPHASE];
    disp(['ph0=' num2str(STARTPHASE)]);
    if STARTPHASE > 1 | STARTPHASE < 0
        warning(['Starting phase = ' num2str(STARTPHASE) ' is outside [0,1]!']);
    end
end

dayInd=[1];
t=[];
lightson=[];
dusk=[];
dawn=[];

dt = 0.1; %must set to <<0.1 because when driving with fractions of 26.7
            %need to be wary of rounding error! May skip dawn or dusk!

t=0:dt:simtime;
for tind=1:numel(t)
    islight = (t(tind) <= dptime+0.5*dt ||...
        t(tind) >= dptime+dpdur+0.5*dt); %dusk is last light pt; dawn is last dark pt
    isdawn = t(tind) >= dptime+dpdur-0.5*dt & t(tind) < dptime+dpdur+0.5*dt;
    isdusk = t(tind) > dptime-0.5*dt & t(tind) <= dptime+0.5*dt;
    disp(['t=' num2str(t(tind)) ...
        '; isdusk=' num2str(isdusk) ...
        '; isdawn=' num2str(isdawn)]);
    
    if islight || tind==1
        ph(tind+1) = ph(tind) + (1/Tlight)*dt;
        lightson(tind) = 1;
    else
        ph(tind+1) = ph(tind) + (1/Tdark)*dt;
        lightson(tind) = 0;
    end
    
    %control always in LL
    phControl(tind+1) = phControl(tind)+(1/Tlight)*dt;
    
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
if (sum(dusk) ~= sum(dawn)) || sum(dusk) ~= 1
    fname = mfilename('fullpath');
    warning([fname ': PO missed a dusk or dawn!']);
end

%%
%find all pts right after lights on, include stepup response
daystart=find(logical(dawn))+1; 
DAWNPHASE = ph(daystart);
if DPTIME + DPDUR > 0+dt & DPTIME + DPDUR < TEND
    DAWNPHASE = ph(logical(dawn));
    DAWNPHASE = mod(DAWNPHASE,1);
    DAWNPHASESHIFT = dawnshift(logical(dawn)); %pshift at last dawn/dusk?
    %DAWNPHASESHIFT = DAWNPHASESHIFT(NUMCYC);
end
 
nightstart=find(logical(dusk))+1;
DUSKPHASE = ph(nightstart);
if DPTIME > 0+dt & DPTIME < TEND
    DUSKPHASE = ph(logical(dusk));
    DUSKPHASE = mod(DUSKPHASE,1);
    DUSKPHASESHIFT = duskshift(logical(dusk));
    %DUSKPHASESHIFT = DUSKPHASESHIFT(NUMCYC);
end

%% plot light-dark boxes overlaid on top of phase oscillator trajectories
if TOMAKELDCYCLES_FIG == 1
    fLDTraj = figure();
    pControl = plot(t,mod(phControl,1),'k.','markersize',6,'linewidth',1);
    hold on;
    
    for n=1:2
        todayLight = lightson & (dayInd == n);
        todayDark = (~lightson) & (dayInd == n);
        pLight(n) = plot(t(todayLight),mod(ph(todayLight),1),'g.','markersize',6,'linewidth',1);
        hold on;
        
        greyCol = [209 210 212]/255;
        yellowCol = [235 235 235]/255;
        
        if sum(todayDark > 0)
            pDark(n) = plot(t(todayDark),mod(ph(todayDark),1),'g.','markersize',6,'linewidth',1);
            pDarkBox = area(t(todayDark),...
                ones(size(t(todayDark))),...
                'facecolor',greyCol,'edgecolor','none');
        end
        
%         %if want to plot this in sine coordinates instead 
%         pLight(n) = plot(t(todayLight),sin(2*pi*(ph(todayLight)-0.5)),'r','linewidth',2);
%         hold on;
%         pDark(n) = plot(t(todayDark),sin(2*pi*(ph(todayDark)-0.5)),'r','linewidth',2);

        if n ~=1
            pLightBox = area(t(todayLight),...
                ones(size(t(todayLight))),...
                'facecolor',yellowCol,'edgecolor','none');
            disp('n~=1');
        else
            disp('n==1');
            pLightBox = area([0 t(todayLight)],...
                [1 ones(size(t(todayLight)))],...
                'facecolor',yellowCol,'edgecolor','none');
        end
        
    end
    
    for n=1:numel(pDark)
        uistack(pDark(n),'top');
    end
    for n=1:numel(pLight)
         uistack(pLight(n),'top');
    end   
    uistack(pControl,'top');
    
%     %picture size [5x3] for illustrator 3 panels
     set(gca,'color','none','layer','top'); %to make sure axes are above colored area bars
     set(gca,'xlim',[0 48],'xtick',0:12:120,'ytick',[0:0.25:1],...
         'xticklabel',[0:12:120],'fontsize',12);
    xlabel('time (hours)');
    ylabel('oscillator phase \theta (rad/2\pi)');
    grid off;
    set(fLDTraj,'units','inches','position',[0 0 4 2.5],'color','w');
    
    if TO_EXPORT_LDCYCLES_FIG == 1
        savename=[getDate('yyyy-mm-dd') '_PhOsc_Illustration_PRC_tDP-' ...
            num2str(DPTIME) '_delta-' num2str(DPDUR) '_' ...
            datestr(now, 'HH.MM.SS') ...
            '.pdf'];
        export_fig(savename,'-cmyk','-painters',fLDTraj);
    end
end

%% plot phase oscillator trajectories, mark dawn and dusk
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
    xlim([0 120]);
    set(gca,'xtick',0:24:120,'ytick',[0:0.25:1]);
    xlabel('hours');
    ylabel('clock phase (rad/2pi)');

end