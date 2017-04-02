%2017-03-28, EL: generate heatmap in Fig. 6E
%simulate entrainment in simple geometric entrainment
%model, assuming infinitely-attracting circular limit cycles with uniform
%speed. 

%dependencies: stepDL.m, stepLD.m, dispif.m
%export figure relies on export_fig.m

%close all;
clear all;
clc;

TODISP=0; %display outputs?
TOEXP=0; %export heatmap?

preXlist=-3:0.3:3;
preRlist=-3:0.3:3;

preRlist = fliplr(preRlist);

Xlist = 10.^preXlist;
Rlist = 10.^preRlist;

OMSIGN = -1; %1 or -1 for CCW vs CW

%save m slopes here
marray=zeros(length(Rlist),length(Xlist));
marray_pk = [];

for x=1:length(Xlist)
    disp(x);
    %disp(['X=' num2str(preXlist(k))]);
    for r=1:length(Rlist)
        dispif(TODISP,['X=' num2str(preXlist(x)) ', R=' num2str(preRlist(r))]);
        X = Xlist(x);
        R = Rlist(r);
        
        T = 24;       
        taulist = 6:1:18;
        entrained = zeros(length(taulist),1);
        
        for i=1:length(taulist)
            tau=taulist(i);
            
            %angle spanned by day vs night
            day = 360 * tau / T;
            night = 360 * (T - tau) / T;
            
            theta = 45;
            ncycles = 30;
            %disp(theta);
            for j = 1:ncycles
                %day
                theta = mod(theta + OMSIGN*day,360);
                %L -> D
                theta = stepLD(theta,R,X);
                %night
                theta = mod(theta + OMSIGN*night,360);
                %D -> L
                theta = stepDL(theta,R,X);
                %disp(theta);
            end

            result = theta; %result= temp stores theta from first simulation
            
            %check independence of starting conditions
            theta = 135;
            ncycles = 30;
            %disp(theta);
            for j = 1:ncycles

                %day
                theta = mod(theta + OMSIGN*day,360);
                %L -> D
                theta = stepLD(theta,R,X);
                %night
                theta = mod(theta + OMSIGN*night,360);
                %D -> L
                theta = stepDL(theta,R,X);
                %disp(theta);
            end
            
            ICdep(r,x) = 1; %by default assume no IC dependence
            
            %check how well two starting conditions agree
            if abs(result - theta) > 1
                marray(r,x) = -0.1;
                marray_pk(r,x) = -0.1;
                ICdep(r,x) = 0;
                dispif(TODISP,'failed IC dependence');
                dispif(TODISP,['reslt=' num2str(result,'%2.2f') ',' ...
                               'theta=' num2str(theta,'%2.2f')]);
            end
            
            entrained_pk(i) = OMSIGN*T*mod((180-theta),360)/360;
            entrained(i) = T*(theta)/360;
        end
        %figure
        %plot(taulist,entrained);
        
        %do linear fit
        coeff = polyfit(taulist', entrained, 1);
        coeff_pk = polyfit(taulist,entrained_pk,1);
        if marray(r,x) ~= -0.1 %save fit only if entrained
            marray(r,x) = coeff(1);
            marray_pk(r,x) = coeff_pk(1);
        end
        
        %get rid of fits that aren't linear
        model = polyval(coeff,taulist');
        fiterror = mean(abs(entrained-model));
        %cutoff
        if fiterror > 0.1*abs(max(entrained)-min(entrained))
            marray(r,x) = -0.1;
            marray_pk(r,x) = -0.1;
            dispif(TODISP,'fails linearity');
            dispif(TODISP,['entr: ' num2str(entrained', '%2.2f ')]);
        end
        dispif(TODISP,['... m=' num2str(marray(r,x))]);
        dispif(TODISP,['... m_pk=' num2str(marray_pk(r,x))]);
        
        %fit better using bestSlope() fn
        posneg = 0;
        [coeff_pk, rsq_pk, breakpt_pk, goodfit_pk] = ...
            bestSlope(taulist', entrained_pk', T, posneg);
        if goodfit_pk == 1 & ICdep(r,x) == 1
            marray_pk(r,x) = coeff_pk;
        end
    end
end

%% make heatmap

fHeatMap=figure();

imagesc(preXlist,preRlist,marray_pk);%,[-0.1 1]);
set(gca,'ydir','normal');
xlabel('log X');
ylabel('log R');
set(gca, 'xtick', -3:1:3,'ytick', -3:1:3);
grid off;

%get current colormap
cmap = colormap;

%make lowest entries = unentrained = white
cmap2 = [1 1 1; cmap(2:end,:)];
colormap(cmap2);
caxis([-0.1 1]);

h=colorbar;
set(h, 'ylim', [0 1]);
xlabel(h,'m');
grid off;

set(fHeatMap,'units','inches','position',[0 0 3.75 3]);
if TOEXP == 1
export_fig([getDate('yyyy-mm-dd') '_mHeatMap_FastRelax_' getDate() '.pdf'],...
    '-cmyk','-painters',fHeatMap);
end
