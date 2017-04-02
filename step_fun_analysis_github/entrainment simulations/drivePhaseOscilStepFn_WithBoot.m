%Drive a phase oscillator with 'light/dark' cycles according to
%stepup/stepdown functions (see code re: response to dawn and dusk). Input
%period of oscillator in light/dark, driving period, fraction of time
%lights are on and the duration of simulation.
%   Input: TDRIVE, TLIGHT, TDARK in hours. DUTYFRAC in [0,1], TEND in hrs.
%   Output: DAWNPHASE [0,1] of last dawn before TEND and corresponding
%   DAWNPHASESHIFT in response to the last dawn. Likewise for DUSK.
%
%
%

function [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
    drivePhaseOscilStepFn_WithBoot(TDRIVE, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
    DUTYFRAC, TEND, NUMCYC, STARTPHASE)
%% TEST
TOTEST=0;
if TOTEST==1
    fname = mfilename('fullpath');
    warning([fname ': Running in test mode!']);
    
    INDIR=['/Users/E/Documents/Advisers/Rust/Photoperiod/'...
        '2017-03-14_mixMatch_LD/analysis'];
    INFILE = '2017-03-17_stepFunMixMatch_inKaiCPhase_19.20.23.mat';
    load([INDIR '/' INFILE]);
    BOOTNUM=1;
    
    TLIGHT=T_hiATP_mix(BOOTNUM);
    TDARK=T_loATP_mix(BOOTNUM);
    STEPUP.phase = up_mix{BOOTNUM}.phase;
    STEPUP.phaseShift = up_mix{BOOTNUM}.phaseShift;
    STEPDOWN.phase = down_mix{BOOTNUM}.phase;
    STEPDOWN.phaseShift = down_mix{BOOTNUM}.phaseShift;
    
    TDRIVE=24.0;
    
    %driving params to both
    DUTYFRAC=[12]/24;
    NUMCYC=1:5;
    TEND=TDRIVE*max(NUMCYC)+1;
    STARTPHASE=0.8;
end

%% make phase oscillator predictions using interpolated stepup/down functions
TODISP=0; %suppress print output during simulations

Tlight = TLIGHT; %upendper (hrs);
Tdark = TDARK; %upstartper (hrs);

Tdrive=TDRIVE;
dutyDrive=DUTYFRAC;
simtime = TEND; %hrs

%get step func. handle
BOOT_STEPFUN = 'getInVitroPhaseShiftBoot';
stepFun = str2func(BOOT_STEPFUN);
dispif(TODISP,['Using ' func2str(stepFun) ' in ' mfilename('fullpath') '.']);

%make figures?
TO_EXPORT_LDCYCLES_FIG = 0;
TOMAKELDCYCLES_FIG = 0;

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

t=double(0:dt:simtime);
Tdrive=double(Tdrive);
tnow=0;

for tind=1:numel(t)
%     tnow = t(tind);
%     while (tnow - Tdrive) > 0
%         tnow = tnow - Tdrive; %mod(tnow,Tdrive);
%     end
%     tnow = mod(tnow,Tdrive);
    
    if tind == 1 || tnow > Tdrive-dt && tnow <= Tdrive+0.5*dt 
        tnow = 0;
    else
        tnow = tnow + dt;
        if tnow > Tdrive-0.5*dt && tnow <= Tdrive+0.5*dt 
            tnow = 0;
        end
    end
    
    islight = tnow <= dutyDrive*Tdrive+0.5*dt &...
        tnow >= 0.5*dt; %dusk is last light pt; dawn is last dark pt
    
    isdawn = tnow >= 0-0.5*dt ... %was >=
        & tnow < 0.5*dt;
    isdusk = tnow > dutyDrive*Tdrive-0.5*dt & ...
        tnow <= dutyDrive*Tdrive+0.5*dt; %was <=
    
%     dispif(TODISP,['t=' num2str(t(tind)) ...
%         '; tnow=' num2str(tnow) ...
%         '; isdusk=' num2str(isdusk) ...
%         '; isdawn=' num2str(isdawn) ...
%         '; islight=' num2str(islight) ...
%         '; Tdrive=' num2str(Tdrive)]);
       
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
        dispif(TODISP, ['dawn at t=' num2str(t(tind)) ',ph=' num2str(ph(tind+1))]);
        %dispif(TODISP, ['dawn pshift=' num2str(stepFun(mod(ph(tind),1),STEPUP))]);
        dawn(tind) = 1;
        dawnshift(tind) = stepFun(mod(ph(tind),1),STEPUP);
        ph(tind+1) = ph(tind+1) + dawnshift(tind);
        dayInd(tind+1) = dayInd(tind)+1;
       
    elseif isdusk & lightson(tind-1) == 1
        dispif(TODISP, ['dusk at t=' num2str(t(tind))]);
        %dispif(TODISP, ['dusk pshift=' num2str(stepFun(mod(ph(tind),1),STEPDOWN))]);
        dusk(tind) = 1;
        duskshift(tind) = stepFun(mod(ph(tind),1),STEPDOWN);
        ph(tind+1) = ph(tind+1) + duskshift(tind);
    end
    
    %detect an error if numdays > numdawns + 1
    if t(tind)/TDRIVE > sum(dawn) + 1 || t(tind)/TDRIVE > sum(dusk) + 1
        error(['Numdays > numdawns or num dusks! Missed a dawn or dusk!']);
    end
    
    %dispif(TODISP,char(10));
end

t(end+1) = t(end)+dt;
lightson=logical(lightson);

dispif(1,['No. dusks: ' num2str(sum(dusk)) ', no. dawns: ' num2str(sum(dawn))]);

%sanity check that haven't skipped a dawn or dusk
if ((sum(dusk) ~= sum(dawn)) || (sum(dusk) ~= max(NUMCYC))) %% && (TODISP == 1)
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
    xlim([0 96]);
    set(gcf,'units','inches','position',[0 0 8.5 3],'color','w');
    
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

