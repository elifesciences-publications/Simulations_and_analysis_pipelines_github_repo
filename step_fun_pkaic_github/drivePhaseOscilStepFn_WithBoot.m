%2017-03-29, EL: drive a phase oscillator with step functions, which are
%passed as arguments. if TOTEST==1, can run this function as a script with
%parameters listed in the if TOTEST==1 {} statement.
%
%Drive a phase oscillator with 'light/dark' cycles according to
%stepup/stepdown functions  (see code re: response to dawn and dusk). Input
%period of oscillator in light/dark, driving period, fraction of time
%lights are on and the duration of simulation.
%   Input: TDRIVE, TLIGHT, TDARK in hours. DUTYFRAC in [0,1], TEND in hrs.
%   Output: DAWNPHASE [0,1] of last dawn before TEND and corresponding
%   DAWNPHASESHIFT in response to the last dawn. Likewise for DUSK.
function [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
    drivePhaseOscilStepFn_WithBoot(TDRIVE, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
    DUTYFRAC, TEND, NUMCYC, STARTPHASE)
%% TEST
TOTEST=0; %run in test mode?
if TOTEST==1
    fname = mfilename('fullpath');
    warning([fname ': Running in test mode!']);
    
    INDIR='saved_data';
    INFILE='2016-10-31_stepFun_Lin-NLin-Oct_1000_2016-10-31_20.45.44.mat';
    load([INDIR '/' INFILE]);
    BOOTNUM=111; %which bootstrapped function to use?
    TLIGHT=T_hiATP(BOOTNUM);
    TDARK=T_loATP(BOOTNUM);
    
    STEPUP.phase = up.phase(BOOTNUM,:);
    STEPUP.phaseShift = up.phaseShift(BOOTNUM,:);
    STEPDOWN.phase = down.phase(BOOTNUM,:);
    STEPDOWN.phaseShift = down.phaseShift(BOOTNUM,:);
    
    TDRIVE=24.0;

    %driving params to both
    DUTYFRAC=[12]/24;
    NUMCYC=1:20;
    TEND=TDRIVE*max(NUMCYC)+1;
    STARTPHASE=0;
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

%explicitly simulate phase oscillator
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
    xlim([0 72]);
    set(gca,'xtick',0:24:120,'ytick',[0:0.25:1]);
    xlabel('hours');
    ylabel('clock phase (rad/2pi)');    
  
end