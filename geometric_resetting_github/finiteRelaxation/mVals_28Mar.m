%2017-03-28, EL: generate Fig 6-figSup2A.
%simulate entrainment in limit cycle geometries for cicular orbits with
%variable attraction strength

%export relies on export_fig.m from FileExchange

close all;
clear all;

INDIR = '.';
cd(INDIR);
OUTDIR = [INDIR '/figs/'];

TOSAVE_ALL = 0; %save workspace variables and timestamped script?

%% simulation and export settings
TOEXP_HEATMAP = 0; %save heatmap?

%select log-space grid of X and R to sample
preX = [-3:0.25:3];
preR = [-3:0.25:3];

preR = fliplr(preR);

%pick values of a to simulate. Note t_half=ln2/a.
aval = 10.^[1 0 -1];
logBase = 10;
Rval = logBase.^preR;
Xval = logBase.^preX;
tauval = [8:2:16]';
NUMCYC = [10]; %how many cycles to entrain to

OMSIGN = 1; %direction of rotation +1 or -1
theta0 = [0 pi];

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
    
        for x=1:numel(Xval)
            disp(['... ... x=' num2str(x)]);
            
            %in manuscript:  
            %day:   rL = 1; xL = 0;
            %night: rD = R; xD = X;
            
            %BUT! in diagrams, havbe
            
            %mimic that 
            dO = makeOrbit(1, 0, aval(a),OMSIGN*2*pi/24);
            nO = makeOrbit(Rval(r), Xval(x), aval(a), OMSIGN*2*pi/24);
            
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
                bestSlopeType='sqerr';
                [MSLOPE_PK(a,r,x,th), RSQ(a,r,x,th), ~, GFIT(a,r,x,th)] = ...
                    bestSlope(tauval,squeeze(PEAKT(a,r,x,th,:)),24,POSNEG,...
                    bestSlopeType);
            end
            for tau = 1:numel(tauval)
                PKDIFF(a,r,x,tau) = ...
                    squeeze(mod(PEAKT(a,r,x,1,tau) - PEAKT(a,r,x,2,tau),24));
                PKDIFF(a,r,x,tau) = wrapVecAround(PKDIFF(a,r,x,tau),12,24,'gt');
            end
            AVGPKDIFF(a,r,x) = mean(abs(squeeze(PKDIFF(a,r,x,:))));
        end
    end
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
    disp(['a=' num2str(a) ', goodpkdiff=']);
    disp(squeeze(goodpkdiff(a,:,:)));
        
    marray_pk(realgood == 0) = -0.1;
    disp(['a=' num2str(a) ', marray_pk=']);
    disp(marray_pk);
    
    %imagesc(preX,preR,marray_pk,[-0.1 1]);
    imagesc(preX,preR,marray_pk);
    title(['t_{1/2}=' num2str(log(2)/aval(a),'%2.2f') ' hr.']); 
    set(gca,'ydir','normal');
    xlabel(['log X']);
    ylabel(['log R']);
    set(gca, 'xtick', -3:1:3,'ytick', -3:1:3);
    grid off;
    
    caxis([-0.1 1]);
    
    %get current colormap
    cmap = colormap(jet);
    cmap2 = [1 1 1; cmap(2:end,:)];
    colormap(cmap2);
    
    h=colorbar;
    set(h, 'ylim', [0 1]);
    xlabel(h,'m');
    grid off;
    
    if TOEXP_HEATMAP == 1
        set(gcf,'units','inches','position',[0 0 2.25 2.25]); % [0 0 4 3.5]);
        export_fig([OUTDIR ...
            getDate('yyyy-mm-dd') '_HeatMap_' 'a' num2str(aval(a)) ...
            '_logBase' num2str(logBase) ...
            '_logX' num2str(min(preX)) '-' num2str(max(preX)) ...
            '_logR' num2str(min(preR)) '-' num2str(max(preR)) ...
            '_' getDate() '.pdf'], ...
            '-cmyk', '-painters', '-pdf', fHeatMap);
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
 
