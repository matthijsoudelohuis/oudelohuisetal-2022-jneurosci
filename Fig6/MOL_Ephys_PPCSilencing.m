%% MOL_PPCSilencing

%% Parameter settings for PSTH
params                      = params_histresponse(); % All time is in microseconds
params.zscore               = 0; %Set z-scoring to zero, normalize baseline later as 100%

params.AlignOn              = 'photostimStart';      %On which timestamp to align as t=0
params.normBaseline         = 1;                %Whether firing rate is normalized to baseline (100%)
params.area                 = 'PPC';
params.clipmodrate          = 2;

%% Get data:
% [Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict' 'ChangeDetectionConflict'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
[Data] = MOL_GetData('E:','CHDET',{'OptoOnly'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
% [Data] = MOL_GetData('E:','CHDET',{'OptoOnly'},[],[],{'sessionData' 'trialData' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Filter sessions:
idx                                 = sessionData.Photostimpower>=2 & strcmp(sessionData.PhotostimArea,'PPC') & strcmp(sessionData.OpsinArea,'PPC');
[sessionData,trialData,spikeData]   = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData,spikeData);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);

fprintf('Subselected %d %s neurons out of %d\n',sum(idx),params.area,length(spikeData.session_ID))
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Main loop to get psth matrix:
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing firing rate for neuron        \n');
snakemat                = NaN(nNeurons,params.nTimebins);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+1));
    fprintf('%d/%d',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %get new trialdata if neuron is recorded in new session:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    hist_mean                       = nanmean(hist_mat,1);
    %Normalize:
    if params.normBaseline
        hist_mean                   = hist_mean/nanmean(hist_mean(params.xtime<0));
    end
    snakemat(iNeuron,:)      = hist_mean;
end
fprintf('\n')

%% Compute modulation ratio:
% modulationratio     = nanmean(snakemat(:,params.xtime>0e6 & params.xtime<1e6),2) ./ nanmean(snakemat(:,params.xtime>-1e6 & params.xtime<-0e6),2);
modulationratio     = nanmean(snakemat(:,params.xtime>0.1e6 & params.xtime<0.9e6),2) ./ nanmean(snakemat(:,params.xtime>-0.9e6 & params.xtime<-0.1e6),2);
modulationratio(modulationratio>params.clipmodrate | isinf(modulationratio)) = params.clipmodrate;

%% Make figure of the distribution over depth:
%Histogram with binning on depth:
edges                   = 0:100:1000;
yticks                  = edges(2:end)-100/2;

%Compute average modulation ratio for inhibited neurons binned on recorded depth:
histY_inh               = NaN(length(yticks),1);
histY_inh_err           = NaN(length(yticks),1);
for iE = 1:length(yticks)-1
    idx                 = spikeData.ChannelY>yticks(iE) & spikeData.ChannelY<=yticks(iE+1) & modulationratio<1;
    histY_inh(iE)       = nanmean(modulationratio(idx));
    histY_inh_err(iE)   = nanstd(modulationratio(idx)) / sqrt(sum(idx));
end

%Compute average modulation ratio for excited neurons binned on recorded depth:
histY_exc               = NaN(length(yticks),1);
histY_exc_err           = NaN(length(yticks),1);
for iE = 1:length(yticks)-1
    idx                 = spikeData.ChannelY>yticks(iE) & spikeData.ChannelY<=yticks(iE+1) & modulationratio>1;
    histY_exc(iE)       = nanmean(modulationratio(idx));
    histY_exc_err(iE)   = nanstd(modulationratio(idx)) / sqrt(sum(idx));
end

%Make figure:
laminardepthfig = figure; set(gcf,'units','normalized','Position',[0.2 0.5 0.24 0.28],'color','w'); hold all;
shadedErrorBar(yticks,histY_inh*100,histY_inh_err*100,{'-b','LineWidth',2},0) %plot inhibited modualtion ratio
shadedErrorBar(yticks,histY_exc*100,histY_exc_err*100,{'-r','LineWidth',2},0) %Plot excited neurons modulation ratio
plot([edges(1) edges(end)],[100 100],'k:','LineWidth',1) %Plot 100% reference line

%Figure make up:
xlabel('Depth from dura (um)','FontSize', 15)
ylabel('Modulation ratio','FontSize',15)
set(gca,'XDir','reverse','linewidth',2)
set(gca,'fontsize',15,'FontWeight','bold','YTick',[50 100 150],'XTick',[0 250 500 750 1000])
xlim([edges(1) edges(end)])
ylim([0 params.clipmodrate*100])
view([90 -90]) %// instead of normal view, which is view([0 90])
text(50,25,sprintf('n=%d inhibited',sum(modulationratio<1)),'FontSize',12)
text(50,110,sprintf('n=%d excited',sum(modulationratio>1)),'FontSize',12)

%% Plot average waveform of excited and inhibited neurons: 
meanwf = cell2mat(spikeData.meanwf');

%Flip neurons that have their spike downward
% Discard neurons that have their spike downward
[nSamples,nNeurons] = size(meanwf);
excludeidx  = abs(min(meanwf(1:30,:))) > abs(max(meanwf(1:30,:)));

%Normalize mean spike waveform to maximum:
for iNeuron = 1:nNeurons
    meanwf(:,iNeuron) = meanwf(:,iNeuron) / max(meanwf(:,iNeuron));
end

%Shift waveforms to align peaks:
params.alignsample      = 16;
for iNeuron = 1:nNeurons
    [~,maxidx]          = max(meanwf(:,iNeuron),[],1);
    sourceidx           = (1:nSamples) + (maxidx-params.alignsample);
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


wv_inh = nanmean(meanwf(:,~excludeidx' & modulationratio<1),2);
wv_exc = nanmean(meanwf(:,~excludeidx' & modulationratio>1),2);

figure; set(gcf,'units','normalized','Position',[0.2 0.5 0.24 0.28],'color','w'); hold all;
plot(wv_inh,'b','LineWidth',2);
plot(wv_exc,'r','LineWidth',2);

figure; set(gcf,'units','normalized','Position',[0.44 0.5 0.24 0.28],'color','w'); hold all;
for iNeuron = 1:nNeurons
    if modulationratio(iNeuron)<1
        plot(meanwf(:,iNeuron),'b','LineWidth',0.5);
    elseif modulationratio(iNeuron)>1
       plot(meanwf(:,iNeuron),'r','LineWidth',0.5);
    end
end


%% 
wv_inh = nanmean(meanwf(:,~excludeidx' & modulationratio<1),2);
wv_exc = nanmean(meanwf(:,~excludeidx' & modulationratio>1),2);

wv_inh_err = nanstd(meanwf(:,~excludeidx' & modulationratio<1),[],2) / sqrt(sum(~excludeidx' & modulationratio<1));
wv_exc_err = nanstd(meanwf(:,~excludeidx' & modulationratio>1),[],2) / sqrt(sum(~excludeidx' & modulationratio>1));

figure; set(gcf,'units','normalized','Position',[0.2 0.5 0.24 0.28],'color','w'); hold all;
xtime = linspace(-0.5,1,48);
shadedErrorBar(xtime,wv_inh,wv_inh_err,{'b','LineWidth',2});
shadedErrorBar(xtime,wv_exc,wv_exc_err,{'r','LineWidth',2});
% xlim([-0.3 1.2])

%Extract features from the average waveform:
for iNeuron = 1:nNeurons
    [iw(iNeuron),ahp(iNeuron),pr(iNeuron),ppd(iNeuron),slope(iNeuron),had(iNeuron)]   =  waveforms_features(meanwf(:,iNeuron)',32000); %#ok<AGROW>
end

% figure; set(gcf,'units','normalized','Position',[0.5 0.5 0.24 0.28],'color','w'); hold all;

meantoplot(1) = nanmean(ppd(~excludeidx' & modulationratio<1));
meantoplot(2) = nanmean(ppd(~excludeidx' & modulationratio>1));
errortoplot(1) = nanstd(ppd(~excludeidx' & modulationratio<1)) / sqrt(sum(~excludeidx' & modulationratio<1));
errortoplot(2) = nanstd(ppd(~excludeidx' & modulationratio>1)) / sqrt(sum(~excludeidx' & modulationratio>1));

% errortoplot(1) = nanstd(ppd(~excludeidx' & modulationratio<1));
% errortoplot(2) = nanstd(ppd(~excludeidx' & modulationratio>1));

% errorbar([1 2],meantoplot,errortoplot)

errorbarxy(meantoplot,[0.45 0.4],errortoplot,errortoplot,[0 0],[0 0])
text(0.7,0.5,'*','FontSize',15)

tempBF = bf.ttest2(ppd(~excludeidx' & modulationratio<1),ppd(~excludeidx' & modulationratio>1));
fprintf('Peak-to-trough delay, %0.2f vs %0.2f ms, %d vs %d neurons, BF=%1.2f\n',meantoplot(1),meantoplot(2),sum(~excludeidx' & modulationratio<1),sum(~excludeidx' & modulationratio>1),tempBF)

% %% Make z-scored figure:
% figure; set(gcf,'units','normalized','Position',[0.1 0.1 0.6 0.5],'color','w')
% % subplot(1,length(splits),iSplit)
% imagesc(snakemat(:,:)*100,params.cscale); hold on;
% plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(snakemat,1)+0.5],'k','LineWidth',5);
% set(gca, 'XTick', 1:1000:length(params.xtime), 'XTickLabels', params.xtime(1:1000:length(params.xtime))/1e6,'FontSize', 20)
% set(gca, 'YTick', [1 size(snakemat,1)], 'YTickLabels', [1 size(snakemat,1)],'FontSize', 20)
% xlabel('Time (s)','FontSize', 20)
% ylabel('Neuron','FontSize', 20)
% c = colorbar;
% colormap(parula);
% c.Label.String = 'Firing rate (% of baseline)';
% 
% %% Make figure of the mean:
% params.colors{1} = [0 0 1];
% 
% figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.6 0.5],'color','w')
% meantoplot = nanmean(snakemat,1) * 100;
% errortoplot = nanstd(snakemat,1)/sqrt(size(snakemat,1)) * 100;
% h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors{1},'LineWidth',3},1);
% h.mainLine.Color = params.colors{1};    h.patch.FaceColor = params.colors{1};
% delete(h.edge(1)); delete(h.edge(2)); 
% handles = h.mainLine; hold all;
% 
% set(gca, 'XTick', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6], 'XTickLabels', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6]./1e6,'FontSize', 20)
% xlim([-0.5e6 1.5e6]);
% ylim([0 120])
% % ylim([0 300])
% 
% ylabel('Firing rate (norm to baseline)','FontSize', 20)
% xlabel('Time (s)','FontSize', 20)

%% 
baseline = nanmean(snakemat(:,params.xtime>-0.9e6 & params.xtime<-0.1e6),2);
resp = nanmean(snakemat(:,params.xtime>0.1e6 & params.xtime<0.9e6),2);

figure; hold all;

scatter(baseline(spikeData.celltype==1),resp(spikeData.celltype==1),'b')
scatter(baseline(spikeData.celltype==2),resp(spikeData.celltype==2),'r')

xlim([0 3])
ylim([0 3])
plot([0 5],[0 5],'k:','LineWidth',1)

% modulationratio(modulationratio>params.clipmodrate | isinf(modulationratio)) = params.clipmodrate;
% 
% scatter(modulationratio

