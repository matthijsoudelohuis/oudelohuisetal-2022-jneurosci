function [IDX_OUT,IDX_LABEL] = calc_neurontype(meanwf)
%IDX_OUT:
%0=unclassified
%1=putative pyramidal
%2=putative interneuron

%%
Par.showFig         = 0;
Par.method          = 'fixedthr'; %{'kmeans' or 'fixedthr'}
% Par.method          = 'kmeans'; %{'kmeans' or 'fixedthr'}
Par.ppd_thres_low   = 0.45; %Peak trough delay lower lim for ppv
Par.ppd_thres_high  = 0.55; %Peak trough delay upper lim for ppyr

colors  = {[0.9 0.1 0.1], [0.1 0.1 0.9], [0.5 0.5 0.5]};
set(0,'Defaultlinelinewidth',2)
set(0,'DefaultAxesLineWidth', 2)

%% Flip neurons that have their spike downward
% Discard neurons that have their spike downward
[~,nNeurons] = size(meanwf);
excludeidx  = abs(min(meanwf(1:30,:))) > abs(max(meanwf(1:30,:)));

IDX_OUT                 = NaN(nNeurons,1);
IDX_OUT(excludeidx)     = 0;

%Remove excluded neurons:
meanwf                  = meanwf(:, ~excludeidx);
[nSamples,nNeurons]     = size(meanwf);

%% Normalize mean spike waveform to maximum:
for iNeuron = 1:nNeurons
    meanwf(:,iNeuron) = meanwf(:,iNeuron) / max(meanwf(:,iNeuron));
end

%% Shift waveforms to align peaks:
params.alignsample      = 16;
for iNeuron = 1:nNeurons
    [~,maxidx]          = max(meanwf(:,iNeuron),[],1);
    sourceidx           = [1:nSamples] + (maxidx-params.alignsample);
    sourceidx           = sourceidx(sourceidx>0 & sourceidx<=nSamples);
    shift               = -(maxidx-params.alignsample);
    
    if shift < 0
        targetidx           = 1:numel(sourceidx);
        meanwf(:,iNeuron)   = [meanwf(sourceidx,iNeuron); interp1(targetidx,meanwf(sourceidx,iNeuron),[1:abs(shift)]+numel(sourceidx),'linear','extrap')'];
    elseif shift > 0
        targetidx           = [1:numel(sourceidx)]+shift;
        meanwf(:,iNeuron)   = [interp1(targetidx,meanwf(sourceidx,iNeuron),[1:abs(shift)],'linear','extrap')'; meanwf(sourceidx,iNeuron)];
    end
    
end

%% Extract Features from the average waveform:
for iNeuron = 1:nNeurons
    [iw(iNeuron),ahp(iNeuron),pr(iNeuron),ppd(iNeuron),slope(iNeuron),had(iNeuron)]   =  waveforms_features(meanwf(:,iNeuron)',32000); %#ok<AGROW>
end

FeatureLabels = {'Initial Width' 'AfterHyperPolarization' 'Peak Ratio' 'Peak-to-peak delay' 'Slope' 'Half-amplitude duration'};
FeatureMat = [iw; ahp; pr; ppd; slope; had;]';

%% sort neurons into groups:
switch Par.method
    case 'kmeans'
        %K-means clustering:
        %Number of clusters: (Pyrs and INs):
        K = 2;
        %Features to use:
        SelectedFeatures = [4 5 6];
        % X = [iw; ppd;had;ahp;pr;slope]'% ahp pr ppd slope had]
        [IDX, clustermeans(:,SelectedFeatures)] = kmeans(FeatureMat(:,SelectedFeatures), K);
        
    case 'fixedthr'
        %Define based solely on peak to trough delay:
        IDX(ppd>Par.ppd_thres_high)=1;
        IDX(ppd<Par.ppd_thres_low)=2;
        IDX(ppd>=Par.ppd_thres_low & ppd<=Par.ppd_thres_high)=3;
end

%% Flip index if necessary to always assign 1 to pyramidals and 2 to INs
if mean(iw(IDX==1)) < mean(iw(IDX==2))
    IDX = 3 - IDX; %Flip IDX
%     clustermeans = flipud(clustermeans);
end

%% Plot figures of the waveform features extracted against each other:
if Par.showFig
    C = nchoosek(1:6,2);
    figure; hold on;
    set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w')
    
    for featurecomb = 1:length(C)
        subplot(4,4,featurecomb); hold on;
        scatter(FeatureMat(IDX==1,C(featurecomb,1)),FeatureMat(IDX==1,C(featurecomb,2)),[],colors{1},'filled'); xlabel(FeatureLabels(C(featurecomb,1)));
        scatter(FeatureMat(IDX==2,C(featurecomb,1)),FeatureMat(IDX==2,C(featurecomb,2)),[],colors{2},'filled'); ylabel(FeatureLabels(C(featurecomb,2)));
        scatter(FeatureMat(IDX==3,C(featurecomb,1)),FeatureMat(IDX==3,C(featurecomb,2)),[],colors{3},'filled');
    end
end

%% Plot figures of the waveform features extracted against each other:

if Par.showFig
    
    %Features to use:
    SelectedFeatures = [4 5 6];
    % X = [iw; ppd;had;ahp;pr;slope]'% ahp pr ppd slope had]
    
    C = nchoosek(SelectedFeatures,2);
    figure; hold on;
    set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w')
    
    for featurecomb = 1:length(C)
        subplot(2,2,featurecomb); hold on;
        scatter(FeatureMat(IDX==1,C(featurecomb,1)),FeatureMat(IDX==1,C(featurecomb,2)),[],colors{1},'filled'); xlabel(FeatureLabels(C(featurecomb,1)));
        %             scatter(clustermeans(1,C(featurecomb,1)),clustermeans(1,C(featurecomb,2)),'k','filled')
        
        scatter(FeatureMat(IDX==2,C(featurecomb,1)),FeatureMat(IDX==2,C(featurecomb,2)),[],colors{2},'filled'); ylabel(FeatureLabels(C(featurecomb,2)));
        %             scatter(clustermeans(2,C(featurecomb,1)),clustermeans(2,C(featurecomb,2)),'k','filled')
        
        scatter(FeatureMat(IDX==3,C(featurecomb,1)),FeatureMat(IDX==3,C(featurecomb,2)),[],colors{3},'filled');
        %     scatter(clustermeans(3,C(featurecomb,1)),clustermeans(3,C(featurecomb,2)),'k','filled')
    end
    
    %         case 'fixedthr'
    %             figure; hold on;
    %             set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w')
    %             edges = [0.1:0.025:1];
    %             bar(edges(1:end-1),histcounts(ppd,edges));
    %
    %     end
    
    
end

%% Show histogram of ppd:
if Par.showFig
    ppd(ppd>=1)=1;
    figure;
    set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w'); hold all;
    binedges = 0:0.025:1.1;
    histY = histcounts(ppd,binedges,'Normalization','count');
    binedges = binedges(1:end-1)+0.0125;
    bar(binedges(binedges>Par.ppd_thres_high),histY(binedges>Par.ppd_thres_high),'r')
    bar(binedges(binedges<Par.ppd_thres_low),histY(binedges<Par.ppd_thres_low),'b')
    bar(binedges(binedges>Par.ppd_thres_low & binedges<Par.ppd_thres_high),histY(binedges>Par.ppd_thres_low & binedges<Par.ppd_thres_high),'k')
    xlim([0.1 1.1]);
    set(gca, 'FontSize', 20);
end


%% Show figure with averaged waveforms and cluster assignment:
if Par.showFig
    figure;
    set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w')
    for iNeuron = 1:nNeurons
        plot(linspace(0,nSamples*0.032,nSamples), meanwf(:,iNeuron),'color', colors{IDX(iNeuron)},'LineWidth',0.5);        hold on
    end
    ylim([-1 1.2]);
    xlim([0 1.3]);
    set(gca, 'FontSize', 20);
    
    figure;
    set(gcf,'units','normalized','Position',[0.2 .2 0.6 0.6],'color','w')
    for i = 1:3
        shadedErrorBar(linspace(0,nSamples*0.032,nSamples),mean(meanwf(:,IDX==i),2),std(meanwf(:,IDX==i),[],2),{'-','LineWidth',1.5,'Color',colors{i}},0); hold on;
    end
    
    ylim([-1 1.2]);
    xlim([0 1.3]);
    set(gca, 'FontSize', 20)
end


%% Verification:

IDX(IDX==3)             = 0; %If any labels are 3 change to 0;
IDX_OUT(~excludeidx)    = IDX;
IDX_LABEL               = repmat({'0=unclassified,1=ppyr,2=ppv'},size(IDX_OUT));

fprintf('Percentage neurons classified as PYR: %2.1f%%\n',sum(IDX==1)/numel(IDX)*100)
fprintf('Percentage neurons classified as PV: %2.1f%%\n',sum(IDX==2)/numel(IDX)*100)
fprintf('Percentage neurons classified as UNCLASSIFIED: %2.1f%%\n',sum(IDX==0)/numel(IDX)*100)

end