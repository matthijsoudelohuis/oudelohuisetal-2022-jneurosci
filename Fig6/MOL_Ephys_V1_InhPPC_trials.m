%% Parameter settings for PSTH
params                      = params_histresponse_coding(); % All time is in microseconds
params.zscore               = 1;
params.conv_win             = 'gaussian';
params.conv_sigma           = 25e3;

%% SnakePlot Parameters:
% params.SortBy           = 'peakLatency';
% params.SortBy           = 'sign';
params.SortBy               = 'maxResponse';
params.AlignOn              = 'stimChange';      %On which timestamp to align as t=0

params.minTrialCond         = 2;

params.area                 = 'V1';
% params.area                 = 'PPC';

%% 
% params.colors_splits       = {[0.2 0.7 0.6] [0.2 0.7 0.6] [0.8 0.2 1] [0.8 0.2 1]};
% params.colors_splits       = {[0.1 0.3 0.5] [0.1 0.3 0.5] [0.4 0.2 0.6]  [0.4 0.2 0.6]};
params.colors_splits       = {[0.1 0.2 1] [0.1 0.2 1] [0.1 0.3 0.5] [0.1 0.3 0.5] [0.8 0 0] [0.8 0 0] [1 0.2 0.4] [1 0.2 0.4]};
params.nSplits              = 8;
params.lines_splits         = {'-' '-.' '-' '-.' '-' '-.' '-' '-.'};
params.labels_splits        = {'VisHit - Ctrl' 'VisHit - Opto' 'VisMiss - Ctrl' 'VisMiss - Opto' 'AuHit - Ctrl' 'AuHit - Opto' 'AuMiss - Ctrl' 'AuMiss - Opto'};

%% Get data:
[Data]              = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
% [Data]              = MOL_GetData('E:','CHDET',{'OptoOnly'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData' 'spikeData'});
% [Data]              = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict'},{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData         = Data.sessionData;
trialData           = Data.trialData;
spikeData           = Data.spikeData;

%% Filter out neurons based on quality:
spikeData           = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter sessions:
idx                                 = sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'PPC') & strcmp(sessionData.OpsinArea,'PPC');
% idx                                 = strcmp(sessionData.PhotostimArea,'PPC') & strcmp(sessionData.OpsinArea,'PPC');
[sessionData,trialData,spikeData]   = MOL_getTempPerSes(sessionData.session_ID(idx),sessionData,trialData,spikeData);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);

for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d sessions, %d trials, %d neurons\n',length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Main loop to compute firing rate:
nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing firing rate for neuron        \n');
snakemat                = NaN(nNeurons,params.nTimebins,params.nSplits);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    splits          = {};
    
    splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==2;
    splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==2;
    
    splits{3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==3;
    splits{4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==3;

    splits{5}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==1;
    splits{6}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==1;
    
    splits{7}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==3;
    splits{8}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==3;
    
    
    
    splits          = {};
    
    splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==2;
    splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==1 & (temptrialData.optoEnd-temptrialData.stimChange)>0.7e6;
    
    splits{3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==3;
    splits{4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.hasphotostim==1 & (temptrialData.optoEnd-temptrialData.stimChange)>0.7e6;

    splits{5}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==1;
    splits{6}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==1;
    
    splits{7}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==0 & temptrialData.vecResponse==3;
    splits{8}       = ismember(temptrialData.audioFreqChangeNorm,[2 3]) & temptrialData.hasphotostim==1 & temptrialData.vecResponse==3;
    
    
    for iSplit = 1:params.nSplits
        if sum(splits{iSplit})>params.minTrialCond
            snakemat(iNeuron,:,iSplit) = nanmean(hist_mat(splits{iSplit},:),1);
        end
    end
    
end

%% Make figure of the mean:
figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.6 0.5],'color','w')

idx_all = true(size(spikeData.session_ID));

handles = [];
for iSplit = 1:4 %params.nSplits
    meantoplot      = nanmean(snakemat(idx_all,:,iSplit),1) * 1;
    errortoplot     = nanstd(snakemat(idx_all,:,iSplit),1)/sqrt(sum(idx_all)) * 1;

    h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit},'markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
    h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
    h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
    handles(iSplit) = h.mainLine; hold all; %#ok<SAGROW>
end

set(gca, 'XTick', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6], 'XTickLabels', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6]./1e6,'FontSize', 20)
xlim([-0.3e6 0.8e6]);

ylabel('Firing rate (normalize to baseline)','FontSize', 20)
xlabel('Time (s)','FontSize', 20)
legend(handles,params.labels_splits); legend boxoff


%% Make figure of the mean:
figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.6 0.5],'color','w')

idx_all = true(size(spikeData.session_ID));

handles = [];
for iSplit = 5:8
    meantoplot      = nanmean(snakemat(idx_all,:,iSplit),1) * 1;
    errortoplot     = nanstd(snakemat(idx_all,:,iSplit),1)/sqrt(sum(idx_all)) * 1;

    if ~all(isnan(meantoplot))
        h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit},'markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
    h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
    h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
    h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
    handles(end+1) = h.mainLine; hold all; %#ok<SAGROW>
    end
end

set(gca, 'XTick', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6], 'XTickLabels', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6]./1e6,'FontSize', 20)
xlim([-0.3e6 0.8e6]);

ylabel('Firing rate (normalize to baseline)','FontSize', 20)
xlabel('Time (s)','FontSize', 20)
legend(handles,params.labels_splits(5:8)); legend boxoff

%% Make figure:
params.colormap = 'parula';

params.cscale = [-0.5 3];

figure; set(gcf,'units','normalized','Position',[0.03 0.45 0.93 0.4],'color','w')
for iSplit = 1:length(splits)
            subplot(1,length(splits),iSplit)
%     subplot(4,2,iSplit)
    imagesc(snakemat(:,:,iSplit),params.cscale); hold on;
    title(params.labels_splits(iSplit));

    plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(snakemat(:,:,iSplit),1)+0.5],'k','LineWidth',5);
    
    set(gca, 'XTick', 1:1000:length(params.xtime), 'XTickLabels', params.xtime(1:1000:length(params.xtime))/1e6,'FontSize', 20)
    set(gca, 'YTick', [1 length(spikeData.ts)], 'YTickLabels', [1 length(spikeData.ts)],'FontSize', 20)
    xlim([find(params.xtime == -0.3e6) find(params.xtime == 1.2e6)]);
    xlabel('Time (s)','FontSize', 20)
%     ylabel('Neuron','FontSize', 20)
    %         c = colorbar;
    colormap(params.colormap);
    c.Label.String = 'Z-scored firing rate';
end

%% Make figure of the hit-miss difference over time: 
snakematdiff                = []; %Construct difference between hits and misses:
snakematdiff(:,:,1)         = diff(snakemat(:,:,[2 1]),[],3); %Diff Max Hit and Miss
snakematdiff(:,:,2)         = diff(snakemat(:,:,[4 3]),[],3);  %Diff Thr Hit and Miss
snakematdiff                = nanmean(snakematdiff,3); %Average over Thr and Max

figure; set(gcf,'units','normalized','Position',[0.3 0.2 0.6 0.5],'color','w')
handles = [];


params.alpha = 0.05/sqrt(params.nTimebins);
params.alpha = 0.05;

idx_UST = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,'VisOnlyTwolevels')));

meantoplot      = nanmean(snakematdiff(idx_UST,:),1) * 1;
errortoplot     = nanstd(snakematdiff(idx_UST,:),1)/sqrt(sum(idx_UST)) * 1;

h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit},'markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
h.mainLine.Color = params.colors_experiments{2};    h.patch.FaceColor = params.colors_experiments{2};
h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
handles(end+1) = h.mainLine; hold all;

statstoplot = NaN(params.nTimebins,1);
ptoplot = NaN(params.nTimebins,1);
for iT = 1:params.nTimebins
    [ptoplot(iT), statstoplot(iT)] = signrank(snakematdiff(idx_UST,iT)',0,'tail','right','alpha',params.alpha);
end
[fdr]        = mafdr(ptoplot); %Correct pvalues for FDR
for iT = 1:params.nTimebins-1
    if fdr(iT)<=params.alpha
        plot([params.xtime(iT) params.xtime(iT+1)],[0.3 0.3],'Color',params.colors_experiments{2},'LineWidth',5)
    end
end

idx_MST = ismember(spikeData.session_ID,sessionData.session_ID(ismember(sessionData.Experiment,'ChangeDetectionConflict')));

meantoplot      = nanmean(snakematdiff(idx_MST,:),1) * 1;
errortoplot     = nanstd(snakematdiff(idx_MST,:),1)/sqrt(sum(idx_MST)) * 1;

h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{params.lines_splits{iSplit},'markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
h.mainLine.Color = params.colors_experiments{3};    h.patch.FaceColor = params.colors_experiments{3};
h.edge(1).Color = [1 1 1];  h.edge(2).Color =  [1 1 1];
h.edge(1).LineWidth = 0.1;  h.edge(2).LineWidth = 0.1;
handles(end+1) = h.mainLine; hold all;

statstoplot = NaN(params.nTimebins,1);
ptoplot = NaN(params.nTimebins,1);
for iT = 1:params.nTimebins
    [ptoplot(iT), statstoplot(iT)] = signrank(snakematdiff(idx_MST,iT)',0,'tail','right','alpha',params.alpha);
end
[fdr]        = mafdr(ptoplot); %Correct pvalues for FDR
for iT = 1:params.nTimebins-1
    if fdr(iT)<=params.alpha
        plot([params.xtime(iT) params.xtime(iT+1)],[0.32 0.32],'Color',params.colors_experiments{3},'LineWidth',5)
    end
end

set(gca, 'XTick', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6], 'XTickLabels', [-0.3e6 0 0.2e6 0.5e6 1e6 1.5e6]./1e6,'FontSize', 20)
xlim([-0.3e6 0.8e6]);
plot(get(gca,'xlim'),[0 0],'k:','LineWidth',2)
ylabel('Firing rate (normalize to baseline)','FontSize', 20)
xlabel('Time (s)','FontSize', 20)
legend(handles,{'Hit-Miss (UST)' 'Hit-Miss (MST)'}); legend boxoff

%% 




