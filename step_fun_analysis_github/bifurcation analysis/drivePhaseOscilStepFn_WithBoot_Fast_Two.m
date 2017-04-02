%2017-03-28, EL: run phase oscillator simulations using experimentally
%measured L and D functions.

%Drive a phase oscillator with 'light/dark' cycles according to
%stepup/stepdown functions. Input
%period of oscillator in light/dark, driving period, fraction of time
%lights are on and the duration of simulation.
%   Input: TDRIVE, TLIGHT, TDARK in hours. DUTYFRAC in [0,1], TEND in hrs.
%   Output: DAWNPHASE [0,1] of last dawn before TEND and corresponding
%   DAWNPHASESHIFT in response to the last dawn. Likewise for DUSK.
%
%Parameter TORUN_FAST toggles between explicit oscillator evolution using
%timetep of 0.01 hrs and calculation that treats L and D as 1-D maps that
%just operate at dawn and dusk.
%
function [DAWNPHASE, DAWNPHASESHIFT, DUSKPHASE, DUSKPHASESHIFT] = ...
    drivePhaseOscilStepFn_WithBoot_Fast_Two(TDRIVE, TLIGHT, TDARK, STEPUP, STEPDOWN, ...
    DUTYFRAC, TEND, NUMCYC, STARTPHASE)
%% TEST
TOTEST=0;
if TOTEST==1
    fname = mfilename('fullpath');
    warning([fname ': Running in test mode!']);
    
    STEPFUNTYPE='BOOT';
    INDIR=['../helper functions and shared files/'];
    INFILE = '2017-04-01_mergedStepFuns_11.35.34.mat';
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

%run fast map iteration or explicit simulation?
TORUN_FAST = 0; %1=fast (preferred), 0=explicit

Tlight = TLIGHT; %upendper (hrs);
Tdark = TDARK; %upstartper (hrs);

Tdrive=TDRIVE;
dutyDrive=DUTYFRAC;
simtime = TEND; %hrs

%get step func. handle
BOOT_STEPFUN = 'getInVitroPhaseShiftBoot';
stepFun = str2func(BOOT_STEPFUN);
dispif(1,['Using ' func2str(stepFun) ' in ' mfilename('fullpath') '.']);

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

%if want to quickly get 1d map iterations
if TORUN_FAST == 1
%go through n cycles
fduskphase=[];
fduskshift=[];
fdawnphase=[];
fdawnshift=[];
tdusk=[];
tdawn=[];

for n=1:numel(NUMCYC)
    if n==1
        fdawnphase(n) = STARTPHASE;
        tdawn(n) = 0;
    else
        fdawnphase(n) = fduskphase(n-1)+fduskshift(n-1)+...
            (1/TDARK)*((1-DUTYFRAC)*TDRIVE);
        fdawnphase(n) = mod(fdawnphase(n),1);
        tdawn(n) = (n-1)*TDRIVE;
    end
    fdawnshift(n) = stepFun(mod(fdawnphase(n),1),STEPUP);
    
    fduskphase(n) = fdawnphase(n) + fdawnshift(n)+...
        (1/TLIGHT)*(DUTYFRAC*TDRIVE);
    fduskphase(n) = mod(fduskphase(n),1);
    tdusk(n) = tdawn(n) + (DUTYFRAC*TDRIVE);
    fduskshift(n) = stepFun(mod(fduskphase(n),1),STEPDOWN);   
    
end
DAWNPHASE = fdawnphase;
DAWNPHASESHIFT = fdawnshift;
DUSKPHASE = fduskphase; 
DUSKPHASESHIFT = fduskshift;
end

%do you want to run explicit simulations with small timestep?
if TORUN_FAST == 0
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
end

%sanity check that haven't skipped a dawn or dusk
if ((numel(DUSKPHASE) ~= numel(DAWNPHASE)) || (numel(DUSKPHASE) ~= max(NUMCYC))) %% && (TODISP == 1)
    fname = mfilename('fullpath');
    warning([fname ': PO missed a dusk or dawn!']);
end

%plot phase oscillator trajectories, mark dawn and dusk
TOPLOT=0;
if TOPLOT==1
    fPhaseOsc = figure();
    plot(t(lightson),mod(ph(lightson),1),'r.','markersize',7,'linewidth',1);
    hold on;
    plot(t(~lightson),mod(ph(~lightson),1),'k.','linewidth',1);
    hold on;
    
    %fdusk, 
    plot(tdusk,mod(fduskphase,1),'c.');
    plot(tdusk,mod(fduskphase+fduskshift,1),'c.');
    plot(tdawn,mod(fdawnphase,1),'y.');
    plot(tdawn,mod(fdawnphase+fdawnshift,1),'y.');
    
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

