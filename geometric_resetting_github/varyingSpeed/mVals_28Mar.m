%2017-03-28, EL: generates Fig. 6-figSup2 B.
%simulations of geometric limit cycle model for orbits with varying angular
%velocity. 
%

close all;
clear all;

INDIR='.';
cd(INDIR);
OUTDIR = [INDIR '/figs/'];

TOSAVE_ALL = 1; %save workspace variables and script with timestamp?

%% simulation and export settings
TOEXP_HEATMAP = 0; %export heatmap?

preX = -3:0.25:3;
preR = -3:0.25:3;

preR = fliplr(preR);

aval = 10.^[0];
preEps =[-6 -2 -1];
epsval = (2.^preEps).*(2*pi/24); %eps. matters relative to omega!

logBase = 10;
Rval = logBase.^preR;
Xval = logBase.^preX;
tauval = [8:2:16]';
NUMCYC = [10]; %how many cycles to entrain to

OMSIGN = 1; %direction of rotation +1 or -1
theta0 = [0 pi];

POSNEG = 0; %for linear fitting; 0 = slope can be both + or -. Else set to +1 or -1.

AVGPKDIFF = [];
MSLOPE_PK = [];
RSQ = [];
PKDIFF = []; %difference in t_pk for different theta0's

%% run simulations
for eps=1:numel(epsval)
    disp(['eps=' num2str(eps)]);
    
    for r=1:numel(Rval)
        disp(['... r=' num2str(r)]);
    
        for x=1:numel(Xval)
            disp(['... ... x=' num2str(x)]);
            
            %in manuscript:  
            %day:   rL = 1; xL = 0;
            %night: rD = R; xD = X;
            
            %BUT! in diagrams, havbe
            
            %mimic that 
            %[orbit] = makeOrbit(R, X, a, eps, omega)
            dO = makeOrbit(1, 0, aval, epsval(eps), OMSIGN*2*pi/24);
            nO = makeOrbit(Rval(r), Xval(x), aval, epsval(eps), OMSIGN*2*pi/24);
            
            for th = 1:numel(theta0)
                disp(['... ... ... th=' num2str(th)]);
                for tau=1:numel(tauval)
                    xyTraj = [];
                    [t, polarTraj, xyTraj, daynight] = ...
                        entrain(dO, nO, tauval(tau), 24, NUMCYC, [dO.R theta0(th)]);
                    [pkT, t2circ, ptOnOrb, appPh0] = ...
                        phaseOnCircle(dO, xyTraj(end,:));
                    
                    PEAKT(eps,r,x,th,tau) = pkT;
                    T2CIRC(eps,r,x,th,tau) = t2circ;
                    PT_ON_ORB(eps,r,x,th,tau,:) = ptOnOrb;
                    APP_PH0(eps,r,x,th,tau) = appPh0;
                    ENDTHETA(eps,r,x,th,tau) = polarTraj(end,2);
                    disp(['tau=' num2str(tau)]);
                end
                
                disp(squeeze(PEAKT(eps,r,x,th,:)));
                bestSlopeType='sqerr';
                [MSLOPE_PK(eps,r,x,th), RSQ(eps,r,x,th), ~, GFIT(eps,r,x,th)] = ...
                    bestSlope(tauval,squeeze(PEAKT(eps,r,x,th,:)),24,...
                    POSNEG,bestSlopeType);
                disp(['goodfit=' num2str(GFIT(eps,r,x,th))]);
            end
            for tau = 1:numel(tauval)
                PKDIFF(eps,r,x,tau) = ...
                    squeeze(mod(PEAKT(eps,r,x,1,tau) - PEAKT(eps,r,x,2,tau),24));
                PKDIFF(eps,r,x,tau) = wrapVecAround(PKDIFF(eps,r,x,tau),12,24,'gt');
            end
            AVGPKDIFF(eps,r,x) = mean(abs(squeeze(PKDIFF(eps,r,x,:))));
        end
    end
end

%% plot heatmaps
for eps=1:numel(epsval)
    fHeatMap=figure();
    
    %tolerate 0.5 hr pk2pk differences if started at different phases
    goodpkdiff = abs(AVGPKDIFF) < 0.5; 
    
    %make sure this isn't transposed relative to fig. in manuscript!
    marray_pk = squeeze(MSLOPE_PK(eps,:,:,1));%MSLOPE(eps,r,x,th)
    good = squeeze(GFIT(eps,:,:,1));
    
    realgood = good & squeeze(goodpkdiff(eps,:,:));
    disp(['eps=' num2str(eps) ', goodpkdiff=']);
    disp(squeeze(goodpkdiff(eps,:,:)));
        
    marray_pk(realgood == 0) = -0.1;
    disp(['eps=' num2str(eps) ', marray_pk=']);
    disp(marray_pk);
    
    %imagesc(preX,preR,marray_pk,[-0.1 1]);
    imagesc(preX,preR,marray_pk);
    title(['\epsilon/\omega =' '2^{' num2str(preEps(eps)) '}']);
    set(gca,'ydir','normal');
    xlabel(['log X']);
    ylabel(['log R']);
    set(gca, 'xtick', -3:1:3,'ytick', -3:1:3);
    grid off;
    
    %get current colormap
    caxis([-0.1 1]);
    cmap = colormap(jet);
    cmap2 = [1 1 1; cmap(2:end,:)];
    colormap(cmap2);
    
%     h=colorbar;
%     set(h, 'ylim', [-0.1 1]);
%     xlabel(h,'m');
    grid off;
    
    if TOEXP_HEATMAP == 1
        set(gcf,'units','inches','position',[0 0 2.25 2.25]);
        export_fig([OUTDIR ...
            getDate('yyyy-mm-dd') '_HeatMap_' 'eps' num2str(epsval(eps),'%2.2f') ...
            '_fitErr-' bestSlopeType ...
            '_logBase' num2str(logBase) ...
            '_logX' num2str(min(preX)) '-' num2str(max(preX)) ...
            '_logR' num2str(min(preR)) '-' num2str(max(preR)) ...
            '_' getDate('HH.MM.SS') '.pdf'], ...
            '-cmyk', '-pdf', fHeatMap);
    end
end

%save workspace and .m file with parameters
if TOSAVE_ALL == 1
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