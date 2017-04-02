%2017-03-28, EL: generate Fig 6-figSup2B.
%simulate entrainment for elliptical limit cycles. 

%dependencies: export_fig.m for figure export, getDate.m to generate
%string timestamps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%These simulations are computationally-intensive and take a 
%long time (~16-20 hrs on a laptop). This code calls functions that 
%explicitly compute the distance between an arbitrary point in space and 
%the nearest point on the elliptical orbit. These numerical solvers, plus
%stiff integration take a long time.
%
% IF YOU ARE MAINLY INTERESTED IN EXAMINING THE RESULTS OF THE SIMULATIONS,
% I HAVE PROVIDED THE OUTPUTS OF THESE SIMULATIONS IN THE /figs/ DIRECTORY 
% ALONG WITH THE SCRIPTS USED TO GENERATE THEM. 
% CONSIDER RUNNING mVals_analyzeSavedWorkspace.m TO LOOK AT THE RESULTS OF
% THE SIMULATIONS.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

INDIR = '.';
OUTDIR = [INDIR '/' 'figs'];
cd(INDIR);

%% simulation and export settings

TOEXP_HEATMAP = 0; %export heatmap?
TOSAVE_ALL = 0; %save all simulation results and the version of the script used to generate them?

%set the ellipticities of the day and night orbits
%eccD = rhoD, eccN = rhoN in Fig. 6-figSup2
eccD = 1; %1 = circle; try 1, 2 10.
eccN = eccD;

%set the range of logX and logR to sample. 
%don't let X==0. start with 0.05 to avoid the situation where one orbit goes
%through the singularity of the other.
preX = 0.05+[-3:0.25:3]; 
preR = [-3:0.25:3];
preR = fliplr(preR);

%set orbit attraction strength
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
                    bestSlope(tauval,squeeze(PEAKT(a,r,x,th,:)),24,POSNEG);
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
        set(gcf,'units','inches','position',[0 0 4 3.5]);
        export_fig([OUTDIR ...
            getDate('yyyy-mm-dd') '_HeatMap_ellipse_' 'a' num2str(aval(a)) ...
            '_logBase' num2str(logBase) ...
            '_logX' num2str(min(preX)) '-' num2str(max(preX)) ...
            '_logR' num2str(min(preR)) '-' num2str(max(preR)) ...
            '_' getDate('HH.MM.SS') '.pdf'], ...
            '-cmyk', '-painters', '-pdf', fHeatMap);
    end
end
    

%save workspace and .m file with parameters
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

