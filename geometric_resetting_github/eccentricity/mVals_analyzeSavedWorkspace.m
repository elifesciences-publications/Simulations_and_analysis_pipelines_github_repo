%2017-03-28, EL: process one set of simulations from ./figs/ folder to plot
%heatmapts for Fig. 6-figSup2B, considering elliptical orbits.
%figure export relies on export_fig.m from FileExchange and getDate.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the following datasets in the workspace to study entrainment
%behaviors for the following elliptical geometries:

%rhoD=1, rhoN=1: 
% mVals_2017-03-15_00.35.36.mat

%rhoD=10, rhoN=10:
% mVals_2017-03-17_00.57.54.mat

%rhoD=10, rhoN=1:
% mVals_2017-03-18_00.17.50.mat

%the code used to generate these datasets is in the mVals_TIMESTAMP.m files
%in the ./figs/ directory, wherE TIMESTAMP corresponds to the timestamp of
%the corresponding .mat file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INDIR = ['.'];
cd(INDIR);
close all;
clear all;

%% simulation and export settings
TOEXP_TRAJ = 0;
TOEXP_HEATMAP = 0;

%% run simulations
for a=1:numel(aval)
    disp(['a=' num2str(a)]);
    
    for r=1:numel(Rval)
        disp(['... r=' num2str(r)]);
    
        for x=1:numel(Xval)
            disp(['... ... x=' num2str(x)]);
                                             
            for th = 1:numel(theta0)
                disp(['... ... ... th=' num2str(th)]);
                                
                disp(squeeze(PEAKT(a,r,x,th,:)));
                [MSLOPE_PK(a,r,x,th), RSQ(a,r,x,th), ~, GFIT(a,r,x,th)] = ...
                    bestSlope(tauval,squeeze(PEAKT(a,r,x,th,:)),24,POSNEG,...
                    'sqerr');
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
    title(['\rho_D=' num2str(eccD) ', \rho_N=' num2str(eccN)]);
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
%     grid off;
    
    TOEXP_HEATMAP = 1;
    if TOEXP_HEATMAP == 1
       
        set(gcf,'units','inches','position',[0 0 2.25 2.25]);
        outfile = [INDIR '/figs/'...
            getDate('yyyy-mm-dd') '_HeatMap_ellipse_' ...
            'a' num2str(aval(a)) ...
            '_eccD' num2str(eccD) ...
            '_eccN' num2str(eccN) ...
            '_logBase' num2str(logBase) ...
            '_logX' num2str(min(preX)) 'to' num2str(max(preX)) ...
            '_logR' num2str(min(preR)) 'to' num2str(max(preR)) ...
            '_' getDate('HH.MM.SS') '.pdf'];
        export_fig(outfile,'-cmyk', '-painters', '-pdf', fHeatMap);
    end
end
    