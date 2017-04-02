%2016-08-23, EL: extract & fit transcriptional dynamics of dawn and dusk
%genes from Vijayan et al., PNAS (2009)

%% load Vikram's processed & annotated data
[vijNum,vijTxt,vijRaw]=xlsread('Vijayan2009_SD1_edited.xlsx','summary');
vijNum = vijNum(2:end,:);

%% load Vikram's raw data (from: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18902)
[rawIntNum,~,~]=xlsread('GSE18902_series_matrix.xlsx','table');

%order raw data in the same way as processed data
order = [];
for n=1:numel(vijNum(:,1))
    if ~isnan(vijNum(n,1))
        order(n,1) = find(rawIntNum(:,1) == vijNum(n,1));
        cleanRawNum(n,:) = rawIntNum(order(n),:);
        cleanVijNum(n,:) = vijNum(n,:);
    end
end

%pick only LL data: GSM468463-79
timepts = [24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84];
cleanRawNum=cleanRawNum(:,1:(1+numel(timepts)));

%%
dawnInd = (cleanVijNum(:,14) == 2);
duskInd = (cleanVijNum(:,14) == 1);
nonCircInd = (cleanVijNum(:,14) == 0);

avgDawn = mean(cleanRawNum(dawnInd,2:end),1);
avgDusk = mean(cleanRawNum(duskInd,2:end),1);
avgNonCirc = mean(cleanRawNum(nonCircInd,2:end),1);

TONORMMINMAX = 0;
if TONORMMINMAX == 1
    normMinMax = @(x) 2*(x-min(x))./(max(x)-min(x)) - 1;
    avgDawn = normMinMax(avgDawn);
    avgDusk = normMinMax(avgDusk);
end

% times for plotting
timesToPlot = 24:0.1:84;

[dawnFit, dawnData] = ...
    fitSinusoidSimple({timepts}, {avgDawn}, 'period',24,'periodLB',23,'periodUB',25);
dawnPred = sinusoidSimple(dawnFit{:}',{timesToPlot},0);

[duskFit, duskData] = ...
    fitSinusoidSimple({timepts}, {avgDusk}, 'period',24,'periodLB',23,'periodUB',25);
duskPred = sinusoidSimple(duskFit{:}',{timesToPlot},0);

fDawnDuskVijayan=figure();
plt(1)=plot(timepts,avgDawn,'bs','markerfacecolor','b');
hold on;
plot(timesToPlot,dawnPred,'b');
plt(2)=plot(timepts,avgDusk,'rs','markerfacecolor','r');
plot(timesToPlot,duskPred, 'r');
plt(3)=plot(timepts,avgNonCirc,'ks','markerfacecolor','k');
set(gca,'xtick',0:6:96,'xlim',[24 84]);
grid off;
legend(plt(1:3),{'dawn genes','dusk genes', 'arrhythmic'});
legend boxoff;

xlabel('time (hours)');
ylabel('norm. expression (AU)');

TOEXPORT_DAWNDUSKVIJAYAN = 0;
if TOEXPORT_DAWNDUSKVIJAYAN == 1
    export_fig('2016-08-23_rawVijayan_Raw_Small_Fit.pdf',...
        '-cmyk','-painters','-pdf',fDawnDuskVijayan);
end

TOSAVEFIT=0;
if TOSAVEFIT==1
    save('2016-08-23_VijayanDawnDuskRawFit.mat','dawnFit','duskFit');
end


%% split up data by Kmeans clusters
pks=[4:4:24];
for km=1:numel(pks)
    kmInd(:,km)=(cleanVijNum(:,13) == pks(km));
    avgKM(km,:)=mean(cleanRawNum(kmInd(:,km),2:end),1);
    errKM(km,:) = std(cleanRawNum(kmInd(:,km),2:end),0,1);
end

fKMeans=figure();
colors = hsv(numel(pks));
for km=1:numel(pks)
    p(km)=plot(timepts,avgKM(km,:),'-','color',colors(km,:));
    erb(km) = errorbar(timepts, avgKM(km,:), errKM(km,:),...
        '-','color',colors(km,:));
    plab{km}=num2str(pks(km));
    hold on;
end
legend(p,plab);
%legend(erb,plab);
set(gca,'xtick',0:4:96);
