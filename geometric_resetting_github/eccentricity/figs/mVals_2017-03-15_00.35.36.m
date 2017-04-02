close all;
clear all;

INDIR = '/home/leypunskiy/eccentricity_Midway_Mar14';
OUTDIR = [INDIR '/' 'figs'];
cd(INDIR);
addpath('/home/leypunskiy/eccentricity_Midway_Mar14/export_fig');
addpath('/home/leypunskiy/MATLAB/');

%% simulation and export settings
TOEXP_TRAJ = 0;
TOEXP_HEATMAP = 1;
TOSAVE_ALL = 1;

%kluge: set up 0.05 + (...) to avoid the situation for X=0, where the one
%cycle goes directly through the singularity of the other; works for
%R=-2:1:2, X=0.05 and theta0=0. will check for X=0.05+(-2:1:2) --> still
%broke for r=1 (=2), x=5 (=2); think this is due to initial condition on
%the X axis that gets confused? tried initial condition th0=(pi/100)+0 and
%exectuted without errors.
%integration takes so long prob. b/c phaseOnCircle always integrates for
%100 hrs! Change to 12 and add in increments of 12 as necessary. --> much
%faster execution time.
%Now try with theta=pi + pi/100.
preX = 0.05+[-3:0.25:3]; %gave error on X=1 for R=-3:1:3, %for X=2, messes  up on R = 3
preR = -3:0.25:3;

%failed on r=6 (preR = -3:0.25:3) and x = 1 (preX = 0.5+[-3:0.25:3]), eccD = eccD =1
%fixed by modifying how circles are detected in ellipseFminbnd -- test if
%ecc = 1/R, as it's defined. not sure why this works better.

%how elliptical are the orbits?
eccD = 1; %1 = circle; try 1, 2 10.
eccN = eccD;

preR = fliplr(preR);

aval = 10.^[1];
logBase = 10;
Rval = logBase.^preR;
Xval = logBase.^preX;
tauval = [8:2:16]';
NUMCYC = [4]; %how many cycles to entrain to

OMSIGN = 1; %direction of rotation +1 or -1
theta0 = (0.01)+[0 pi];

POSNEG = 0; %for linear fitting; 0 = slope can be + or -. Else set to +1 or -1.

AVGPKDIFF = [];
MSLOPE_PK = [];
RSQ = [];
PKDIFF = []; %difference in t_pk for different theta0's

%% run simulations
for a=1:numel(aval)
    disp(['a=' num2str(a)]);
    
    for r=1:numel(Rval)
        disp(['... r=' num2str(r)]);
        
        fTestOrb = [];
        for x=1:numel(Xval)
            disp(['... ... x=' num2str(x)]);
            
            %in manuscript:  
            %day:   rL = 1; xL = 0;
            %night: rD = R; xD = X;
            
            %BUT! in diagrams, havbe
            %[orbit] = makeEllipse(R, X, delta, attraction, omega)
            %must have R*delta < ecc. ecc = 1 --> circle.
            
            %mimic that 
            dO = makeEllipse(1, 0, 1/eccD, aval(a),OMSIGN*2*pi/24);
            nO = makeEllipse(Rval(r), Xval(x), 1/(eccN*Rval(r)), aval(a), OMSIGN*2*pi/24);
            
            % plot orbits
            close(fTestOrb);
            fTestOrb = figure();
            [~,~,dCircle] = ellipseTraj(dO);
            [~,~,nCircle] = ellipseTraj(nO);
            plot(dCircle(:,1),dCircle(:,2),'r','linewidth',1);
            hold on;
            plot(nCircle(:,1),nCircle(:,2),'k','linewidth',1);
            title(['R=' num2str(logBase)  '^' num2str(preR(r)) ', ' ...
                'X=' num2str(logBase) '^' num2str(preX(x))]);
            axis('equal');
                      
            for th = 1:numel(theta0)
                disp(['... ... ... th=' num2str(th)]);
                for tau=1:numel(tauval)
                    xyTraj = [];
                    [t, polarTraj, xyTraj, daynight] = ...
                        entrain(dO, nO, tauval(tau), 24, NUMCYC, [dO.R theta0(th)]);
                    [pkT, t2circ, ptOnOrb, appPh0] = ...
                        phaseOnCircle(dO, xyTraj(end,:));
                    
                    PEAKT(a,r,x,th,tau) = pkT;
                    T2CIRC(a,r,x,th,tau) = t2circ;
                    PT_ON_ORB(a,r,x,th,tau,:) = ptOnOrb;
                    APP_PH0(a,r,x,th,tau) = appPh0;
                    ENDTHETA(a,r,x,th,tau) = polarTraj(end,2);
                end
                
                disp(squeeze(PEAKT(a,r,x,th,:)));
                [MSLOPE_PK(a,r,x,th), RSQ(a,r,x,th), ~, GFIT(a,r,x,th)] = ...
                    bestSlope(tauval,squeeze(PEAKT(a,r,x,th,:)),24,POSNEG,'sqerr');
            end
            if numel(theta0) > 1
                for tau = 1:numel(tauval)
                    PKDIFF(a,r,x,tau) = ...
                        squeeze(mod(PEAKT(a,r,x,1,tau) - PEAKT(a,r,x,2,tau),24));
                    PKDIFF(a,r,x,tau) = wrapVecAround(PKDIFF(a,r,x,tau),12,24,'gt');
                end
                AVGPKDIFF(a,r,x) = mean(abs(squeeze(PKDIFF(a,r,x,:))));
            end
        end
    end
end

%% save workspace and .m file with parameters
if TOSAVE_ALL
    %get current .m file name
    FileNameAndLocation=[mfilename('fullpath')];
    NameOnly = strsplit(FileNameAndLocation,'/');
    NameOnly = NameOnly{end};
    
    %append date to file name
    NewName = [NameOnly '_' getDate()]; 
    currentfile=strcat(FileNameAndLocation, '.m');
    
    %copy file with new name
    copyfile(currentfile,[OUTDIR '/' NewName '.m']);
    
    %save workspace with the same name
    save([OUTDIR '/' NewName '.mat'])
end

%% plot heatmaps
for a=1:numel(aval)
    fHeatMap=figure();
    
    %tolerate 0.5 hr pk2pk differences if started at different phases
    goodpkdiff = abs(AVGPKDIFF) < 0.5; 
    
    %make sure this isn't transposed relative to fig. in manuscript!
    
    marray_pk = squeeze(MSLOPE_PK(a,:,:,1));%MSLOPE(a,r,x,th)
    good = squeeze(GFIT(a,:,:,1));
    
    realgood = good & squeeze(goodpkdiff(a,:,:));
    
    marray_pk(realgood == 0) = -0.1;
    
    %imagesc(preX,preR,marray_pk,[-0.1 1]);
    imagesc(preX,preR,marray_pk);
    title(['R/\delta_D=' num2str(eccD) ', ' 'R/\delta_N=' num2str(eccN)]);
    set(gca,'ydir','normal');
    xlabel(['log X']); % (base ' num2str(logBase) ')']);
    ylabel(['log R']); % (base ' num2str(logBase) ')']);
    set(gca, 'xtick', -3:1:3,'ytick', -3:1:3);
    grid off;
    
    %same color axis scale
    caxis([-0.1 1]);
    
    %get current colormap
    cmap = colormap;
    cmap2 = [1 1 1; cmap(2:end,:)];
    colormap(cmap2);
    
    h=colorbar;
    set(h, 'ylim', [-0.1 1]);
    xlabel(h,'m');
    grid off;
    
    if TOEXP_HEATMAP == 1
        set(gcf,'units','inches','position',[0 0 2.25 2.25]);
        export_fig([OUTDIR '/'...
            getDate('yyyy-mm-dd') '_HeatMap_ellipse_' 'a' num2str(aval(a)) ...
            '_logBase' num2str(logBase) ...
            '_logX' num2str(min(preX)) '-' num2str(max(preX)) ...
            '_logR' num2str(min(preR)) '-' num2str(max(preR)) ...
            '_' getDate('HH.MM.SS') '.pdf'], ...
            '-cmyk', '-painters', '-pdf', fHeatMap);
    end
end 

break;

%% plot trajectories
nightCol = [0.2 0.2 0.2];
dayCol = [255 215 0]./255;
pickA = 1;
pickR = 4;
pickX = 4;
pickTheta = 2;
%pickTau = 5;

SimVal.aval = aval;
SimVal.Rval = Rval;
SimVal.Xval = Xval;
SimVal.tauval = tauval;
SimVal.OMSIGN = OMSIGN;
SimVal.NUMCYC = NUMCYC;
SimVal.theta0 = theta0;
SimVal.PEAKT = PEAKT;
SimVal.DAY_DELTA = 1/2;
SimVal.NIGHT_DELTA = 1/(2*Rval(pickR));

SimPick.pickA = pickA;
SimPick.pickR = pickR;
SimPick.pickX = pickX;
%SimPick.pickTau = pickTau;
SimPick.pickTheta = pickTheta;

%% plot static trajectory
fTraj = figure('color',[0.9 0.9 0.9]);
cnt = 1;
pickTau = [1 3 5];
for p = 1:numel(pickTau)
    SimPick.pickTau = pickTau(p);
    subplot(2*numel(pickTau)+1,1,(2*cnt):(2*cnt+1));
    [a,b,c] = plotEntrain(SimPick,SimVal,24*[NUMCYC-1 NUMCYC],1);
    cnt=cnt+1;
end

figTit = ...
    ['a=' num2str(aval(pickA)) ...
    '. logX=' num2str(preX(pickX)) ...
    ', logR=' num2str(preR(pickR)) ...
    ' (base ' num2str(logBase) ')' 10 ...
    '2\pi/\omega=' num2str(2*pi/dO.omega) ' hr' ...
    ', \theta_0=' num2str(theta0(pickTheta))];
suplabel(figTit,'t'); %,'t',[.08 .07 .84 .84]);
set(gcf,'units','inches','position',[0 0 5 6]);


TOEXP_TRAJ = 1;
if TOEXP_TRAJ == 1
    OUTDIR = ['/Users/E/Documents/Advisers/Rust/Simulations/' ... 
        'geometric resetting/eccentricity/figs/'];
    export_fig([OUTDIR ...
        getDate('yyyy-mm-dd') '_xyTraj_' 'a' num2str(aval(pickA)) ...
        'logBase' num2str(logBase) '_logX' num2str(preX(pickX)) ...
        '_logR'  num2str(preR(pickR)) ...
        '_th' num2str(theta0(pickTheta)) '_' getDate('HH.MM.SS') '.pdf'], ...
        '-cmyk', '-painters', '-pdf', fTraj);
end
