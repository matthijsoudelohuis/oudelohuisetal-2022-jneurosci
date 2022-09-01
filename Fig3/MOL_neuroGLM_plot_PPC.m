%% MOL_neuroGLM_plot

%% Load dataset:
% savedate                    = '09-07-20';
% savedate                    = '04-01-21';
% savedate                    = '14-12-20';
savedate                    = '19-03-21';
folderpath                  = fullfile('E:','Data','Analysis','neuroGLM',savedate);
load(fullfile(folderpath,'neuroGLM.mat'))
params.savedir              = 'E:\Documents\PhD\Figures\Project CHDET\Results - PPC\Figure 2 - GLM mixed encoding';
params.colorsplits{5}       = [0 0.55 0.2];

%% figure settings:
set(0,'DefaultAxesLineWidth',1)
set(0,'defaultLineLineWidth',1)
set(0,'defaultaxesfontsize',10)

%% Output figure parameters
params.nMaxTrials           = max(trialData.trialNum);
params.nPredictors          = size(output.x,3);

params.Experiments          = {'ChangeDetectionConflictDecor' 'VisOnlyTwolevels' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels     = {'NE' 'UST' 'MST'}; %Which versions of the task to load data from
params.ExperimentLabels     = {'Naive' 'XXXX' 'Trained'}; %Which versions of the task to load data from

params                      = MOL_getColors_CHDET(params);

%% Compute variance explained over all individual trials:
idx         = params.xtime>0 & params.xtime<0.5e6;
idx         = repmat(idx,1,params.nMaxTrials);
for iNeuron = 1:params.nNeurons %loop over neurons
    for iM = 1:params.nModels %loop over models
        %Compute variance explained by the model for this neuron (across certain time range)
        output.var_expl_v1(iNeuron,iM)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat(iNeuron,iM,idx))) / nanvar(output.y(iNeuron,idx));
%         output.var_expl_v1(iNeuron,iM)     = nanvar(squeeze(output.y_hat(iNeuron,iM,idx))) / nanvar(output.y(iNeuron,idx));
    end
end

%% Compute variance explained over trialtype-averaged activity (Runyan et al. Harvey lab custom)
idx         = params.xtime>0 & params.xtime<0.5e6;
for iNeuron = 1:params.nNeurons %loop over neurons
    [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
    nTotalSpikeBins     = sum(~isnan(output.y(iNeuron,:)));
    nTrials             = nTotalSpikeBins/params.nTimebins;
    hist_mat            = reshape(output.y(iNeuron,1:nTotalSpikeBins),params.nTimebins,nTrials);
    
    splitidxs{1} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{2} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
    splitidxs{3} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==3;
    splitidxs{4} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==1;
%     splitidxs{5} = strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==3;
    
    hist_mat_mean = [];
    for iSplit = 1:length(splitidxs)
        hist_mat_mean = [hist_mat_mean; nanmean(hist_mat(idx,splitidxs{iSplit}),2)]; %#ok<*AGROW>
    end
    
    for iM = 1:params.nModels
        hist_model          = reshape(output.y_hat(iNeuron,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
        
        hist_model_mean = [];
        for iSplit = 1:length(splitidxs)
            hist_model_mean = [hist_model_mean; nanmean(hist_model(idx,splitidxs{iSplit}),2)];
        end
        output.var_expl_v2(iNeuron,iM)     = 1 - nanvar(hist_mat_mean - hist_model_mean) / nanvar(hist_mat_mean);
    end
end

%% Figure variance explained:
figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.05 0.4 0.25],'color','w');
for iExp = 1:length(params.Experiments)
    subplot(1,length(params.Experiments),iExp); hold all;
    title(params.ExperimentLabels(iExp))
    expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    violinplot(output.var_expl_v1(expidx,:),[],'Width',0.45,'ViolinColor',params.colors_experiments{iExp},'ViolinAlpha',1);
    set(gca,'xticklabels',params.modelString)
    ylim([0 max(output.var_expl_v1(:))]);
    ylabel('Explained variance')
end

%% Figure variance explained:
%Figure variance explained:
figure; hold all; set(gcf,'units','normalized','Position',[0.45 0.35 0.16 0.32],'color','w');
exp_selec = [1 3];
for iExp = 1:length(exp_selec) %length(params.Experiments)
    expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(exp_selec(iExp)))));
    errorbar(iExp,nanmean(output.var_expl_v2(expidx,2)),nanstd(output.var_expl_v2(expidx,2)),'k.','LineWidth',2);
    h = bar(iExp,nanmean(output.var_expl_v2(expidx,2)),0.45,'k');
    set(h,'FaceColor',params.colors_experiments{exp_selec(iExp)});
end
set(gca,'xtick',[1 2],'xticklabels',params.ExperimentLabels(exp_selec))
ylim([0 1])
ylabel('Explained variance')

iExp = 1;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
xdata = output.var_expl_v2(expidx,2);
iExp = 3;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
ydata = output.var_expl_v2(expidx,2);

tempBF = bf.ttest2(xdata,ydata);
bayesstar([1 2],tempBF);

%% 
fprintf('Percent variance explained trial-type average (MST): %2.1f%%\n',nanmean(output.var_expl_v2(expidx,2)*100))
fprintf('Percent variance explained trial to trial (MST): %2.1f%%\n',nanmean(output.var_expl_v1(expidx,2)*100))

%% Compute variance explained over all individual trials for subselection of variables:
idx         = params.xtime>0 & params.xtime<0.5e6;
idx         = repmat(idx,1,params.nMaxTrials);

var_expl_splits       = NaN(params.nNeurons,params.nSplits);
for iS = 1:params.nSplits
    for iNeuron = 1:params.nNeurons %loop over neurons
        %Compute variance explained by this split (subselection of predictors) for this neuron (across certain time range)
        var_expl_splits(iNeuron,iS)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat_split(iNeuron,iS,idx))) / nanvar(output.y(iNeuron,idx));
%         var_expl_splits(iNeuron,iS)     = nanvar(squeeze(output.y_hat_split(iNeuron,iS,idx))) / nanvar(output.y(iNeuron,idx));
    end
end

%Exclude infinite value occuring in a neuron;
var_expl_splits(var_expl_splits>1 | var_expl_splits<-1) = NaN;

%% Create figure PSTH:
% temp                        = output.var_expl_v1(:,2) - output.var_expl_v1(:,1);
% temp                        = output.var_expl_v2(:,2) - output.var_expl_v2(:,1);
% params.example_cell_IDs     = []; %init var to store some cell ids of example neurons
% for iExp = 3%1:length(params.Experiments) %take some examples from each experiment type:
%     expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
%     params.example_cell_IDs = [params.example_cell_IDs; spikeData.cell_ID(temp > prctile(temp(expidx),98,1) & expidx)]; %with high expl variance
% end

iM = 2; %show only for full model:

h = spikeData.cell_ID(expidx);
% params.example_cell_IDs = h(sigvarmat(1,:) & ~sigvarmat(2,:) & ~sigvarmat(3,:));

params.example_cell_IDs = {
    '20102018081321110'
    '20112018081021033'
    '20092018082321319'
    '20302020011611111'
    '20122018081421187'
    '20122018081421163'
    '20312020011511420'
    };


% params.example_cell_IDs = spikeData.cell_ID(1:10);
colors      = {[12 168 237]      [96 12 237]         [255 140 0]     [232 28 12]};
colors      = cellfun(@(x) x/256,colors,'UniformOutput',false);
labels      = {'Visual Miss'    'Visual Hit'    'Audio Miss'    'Audio Hit'};

nExNeurons = length(params.example_cell_IDs);

figure; hold all;
set(gcf,'units','normalized','Position',[0.1 0.08 0.3 nExNeurons*0.136],'color','w')

for iNeuron = 1:nExNeurons
    
    neuron_idx = strcmp(spikeData.cell_ID,params.example_cell_IDs(iNeuron));
    if any(neuron_idx)
        [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(neuron_idx),trialData);
        
        nTotalSpikeBins     = sum(~isnan(output.y(neuron_idx,:)));
        nTrials             = nTotalSpikeBins/params.nTimebins;
        hist_mat            = reshape(output.y(neuron_idx,1:nTotalSpikeBins),params.nTimebins,nTrials);
        hist_model          = reshape(output.y_hat(neuron_idx,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
        
        %Show for these trial types the avg firing rate and model estimated firing rate:
        splitidxs = {};
        splitidxs{1} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
        splitidxs{2} = strcmp(temptrialData.trialType,'X') & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
        splitidxs{3} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==3;
        splitidxs{4} = strcmp(temptrialData.trialType,'Y') & temptrialData.audioFreqChangeNorm==3 & temptrialData.vecResponse==1;
        %     splitidxs{5} = strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==3;
        
        nSplits = length(splitidxs);
        
        hist_mat_mean           = NaN(nSplits,params.nTimebins);
        hist_mat_model_mean     = NaN(nSplits,params.nTimebins);
        
        for iSplit = 1:nSplits
            %Get the mean over trials::
            hist_mat_mean(iSplit,:)           = nanmean(hist_mat(:,splitidxs{iSplit}),2);
            hist_mat_model_mean(iSplit,:)      = nanmean(hist_model(:,splitidxs{iSplit}),2);
        end
        
        subplot(nExNeurons,2,(iNeuron-1)*2+1); hold all;
        CurrLineHandles = [];
        for iSplit = 1:2
            
            CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean(iSplit,:),'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
            plot(params.xtime,hist_mat_model_mean(iSplit,:),':','LineWidth',2,'Color',colors{iSplit});
            xlim([params.xtime(1) params.xtime(end)]);
            ylabel(''); xlabel('');
            set(gca,'XTick',[],'YTick',[]);
            plot([0 0],[0 100],'k:','LineWidth',2)
            
            ylim([0 ceil(max(hist_mat_mean(:)))])
            set(gca,'YTick',ceil(max(hist_mat_mean(:))))
        end
        RT = nanmean(temptrialData.responseLatency(splitidxs{2}));
        plot(RT,ceil(max(hist_mat_mean(:)))*0.95,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
        
        subplot(nExNeurons,2,(iNeuron-1)*2+2); hold all;
        for iSplit = 3:4
            
            CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean(iSplit,:),'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
            plot(params.xtime,hist_mat_model_mean(iSplit,:),':','LineWidth',2,'Color',colors{iSplit});
            
            xlim([params.xtime(1) params.xtime(end)]);
            ylabel(''); xlabel('');
            set(gca,'XTick',[],'YTick',[]);
            plot([0 0],[0 100],'k:','LineWidth',2)
            ylim([0 ceil(max(hist_mat_mean(:)))])
            set(gca,'YTick',ceil(max(hist_mat_mean(:))))
            
        end
        
        RT = nanmean(temptrialData.responseLatency(splitidxs{4}));
        plot(RT,ceil(max(hist_mat_mean(:)))*0.95,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
        
        %Figure makeup:
        %         xlabel('Time (sec)','FontSize',15)    %Label x-axis
        %     ylabel('Hz','FontSize',15)            %Label y-axis
        %         set(gca, 'XTick', params.xtime(1:500e3/params.binsize:params.nTimebins), 'XTickLabels', params.xtime(1:500e3/params.binsize:params.nTimebins)/1e6,'FontSize', 20)
        %     set(gca, 'FontSize', 20);
        %     legend(CurrLineHandles,labels,'FontSize',15);
        %     legend(gca,'boxoff');
        %     end
    end
end
tightfig(gcf);

%% 
% iExp = 3;
% expidx          = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
% 
% varidx          = ~cellfun(@isempty,strfind(params.varsplits,'history'));
% thresh          = prctile(var_expl_splits(expidx,varidx),95);
% 
% params.example_cell_IDs = spikeData.cell_ID(var_expl_splits(:,varidx)>thresh & expidx);

%% Figure example cell history effect:
params.example_cell_IDs = {'20092018082321322'};

colors      = {[0.5 0 0 ]       [1 0.4 0.05]};
labels      = {'No Reward'         'Reward'};


nExNeurons = length(params.example_cell_IDs);

figure; hold all;
set(gcf,'units','normalized','Position',[0.4 0.08 0.3 nExNeurons*0.136],'color','w')

for iNeuron = 1:nExNeurons
    
    neuron_idx = strcmp(spikeData.cell_ID,params.example_cell_IDs(iNeuron));
    if any(neuron_idx)
    [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(neuron_idx),trialData);
    
    nTotalSpikeBins     = sum(~isnan(output.y(neuron_idx,:)));
    nTrials             = nTotalSpikeBins/params.nTimebins;
    hist_mat            = reshape(output.y(neuron_idx,1:nTotalSpikeBins),params.nTimebins,nTrials);
    hist_model          = reshape(output.y_hat(neuron_idx,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
    
    %Show for these trial types the avg firing rate and model estimated firing rate:
    splitidxs = {};
    splitidxs{1} = [false; cellfun(@isempty,temptrialData.rewardTime(1:end-1))];% & temptrialData.vecResponse(1:end-1)~=3];
    splitidxs{2} = ~splitidxs{1};
    
%     splitidxs{1} = [false; temptrialData.correctResponse(1:end-1)==1 & temptrialData.vecResponse(1:end-1)~=3];
%     splitidxs{2} = ~splitidxs{1};
                
    nSplits = length(splitidxs);
    
    hist_mat_mean           = NaN(nSplits,params.nTimebins);
    hist_mat_model_mean     = NaN(nSplits,params.nTimebins);
    
    for iSplit = 1:nSplits
        %Get the mean over trials::
        hist_mat_mean(iSplit,:)           = nanmean(hist_mat(:,splitidxs{iSplit}),2);
        hist_mat_model_mean(iSplit,:)      = nanmean(hist_model(:,splitidxs{iSplit}),2);
    end
    
    
    iM = 2; %show only for full model:
        
    subplot(nExNeurons,2,(iNeuron-1)*2+1); hold all;
    CurrLineHandles = [];
    for iSplit = 1:2
        
        CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean(iSplit,:),'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
        plot(params.xtime,hist_mat_model_mean(iSplit,:),':','LineWidth',2,'Color',colors{iSplit});
        xlim([params.xtime(1) params.xtime(end)]);
        ylabel(''); xlabel('');
        set(gca,'XTick',[],'YTick',[]);
        plot([0 0],[0 100],'k:','LineWidth',2)
        
        ylim([0 ceil(max(hist_mat_mean(:)))])
        set(gca,'YTick',ceil(max(hist_mat_mean(:))))
    end
        %     Figure makeup:
%         xlabel('Time (sec)','FontSize',15)    %Label x-axis
%         ylabel('Hz','FontSize',15)            %Label y-axis
%         set(gca, 'XTick', params.xtime(1:500e3/params.binsize:params.nTimebins), 'XTickLabels', params.xtime(1:500e3/params.binsize:params.nTimebins)/1e6,'FontSize', 20)
        set(gca, 'FontSize', 20);
        legend(CurrLineHandles,labels,'FontSize',15);
        legend(gca,'boxoff');
    end
end
tightfig(gcf);

%% Figure example cell history effect:
params.example_cell_IDs = {'20092018082321324'};

colors      = {[1 0.4 0.05]      [0.5 0 0]};
labels      = {'Previous Choice Au'   'Previous Choice Vis'};

nExNeurons = length(params.example_cell_IDs);

figure; hold all;
set(gcf,'units','normalized','Position',[0.4 0.08 0.3 nExNeurons*0.136],'color','w')

for iNeuron = 1:nExNeurons
    
    neuron_idx = strcmp(spikeData.cell_ID,params.example_cell_IDs(iNeuron));
    if any(neuron_idx)
        [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(neuron_idx),trialData);
        
        nTotalSpikeBins     = sum(~isnan(output.y(neuron_idx,:)));
        nTrials             = nTotalSpikeBins/params.nTimebins;
        hist_mat            = reshape(output.y(neuron_idx,1:nTotalSpikeBins),params.nTimebins,nTrials);
        hist_model          = reshape(output.y_hat(neuron_idx,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
        
        %Show for these trial types the avg firing rate and model estimated firing rate:
        splitidxs = {};
        
        splitidxs{1} = [false; temptrialData.vecResponse(1:end-1)==1];
        splitidxs{2} = [false; temptrialData.vecResponse(1:end-1)==2];
        
        nSplits = length(splitidxs);
        
        hist_mat_mean           = NaN(nSplits,params.nTimebins);
        hist_mat_model_mean     = NaN(nSplits,params.nTimebins);
        
        for iSplit = 1:nSplits
            %Get the mean over trials::
            hist_mat_mean(iSplit,:)           = nanmean(hist_mat(:,splitidxs{iSplit}),2);
            hist_mat_model_mean(iSplit,:)      = nanmean(hist_model(:,splitidxs{iSplit}),2);
        end
        
        
        iM = 2; %show only for full model:
        
        subplot(nExNeurons,2,(iNeuron-1)*2+1); hold all;
        CurrLineHandles = [];
        for iSplit = 1:2
            
            CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean(iSplit,:),'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
            plot(params.xtime,hist_mat_model_mean(iSplit,:),':','LineWidth',2,'Color',colors{iSplit});
            xlim([params.xtime(1) params.xtime(end)]);
            ylabel(''); xlabel('');
            set(gca,'XTick',[],'YTick',[]);
            plot([0 0],[0 100],'k:','LineWidth',2)
            
            ylim([0 ceil(max(hist_mat_mean(:)))])
            set(gca,'YTick',ceil(max(hist_mat_mean(:))))
        end
        %     Figure makeup:
        %         xlabel('Time (sec)','FontSize',15)    %Label x-axis
        %         ylabel('Hz','FontSize',15)            %Label y-axis
                set(gca, 'XTick', [0 1.5e6], 'XTickLabels', [0 1.5e6]/1e6,'FontSize', 20)
        set(gca, 'FontSize', 20);
        legend(CurrLineHandles,labels,'FontSize',15);
        legend(gca,'boxoff');
    end
end
tightfig(gcf);

%% Figure example cell lick modulation:
iExp = 3;
expidx          = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

varidx          = ~cellfun(@isempty,strfind(params.varsplits,'motor'));
thresh          = prctile(var_expl_splits(expidx,varidx),80);

% params.example_cell_IDs = spikeData.cell_ID(var_expl_splits(:,varidx)>thresh & expidx);

params.example_cell_IDs = {'20352019121921873'};

colors      = {[0.5 0 0]    [1 0.4 0.05]};
labels      = {'Lick'   'No Lick'};

nExNeurons = length(params.example_cell_IDs);

figure; hold all;
set(gcf,'units','normalized','Position',[0.4 0.08 0.3 nExNeurons*0.136],'color','w')
% set(gcf,'units','normalized','Position',[0.4 0.08 0.3 1],'color','w')

% params.example_cell_IDs([18 19 20])
% params.example_cell_IDs([4 9])
for iNeuron = 1:nExNeurons
    
    neuron_idx = strcmp(spikeData.cell_ID,params.example_cell_IDs(iNeuron));
    if any(neuron_idx)
    [temptrialData]     = MOL_getTempPerSes(spikeData.session_ID(neuron_idx),trialData);
    
    nTotalSpikeBins     = sum(~isnan(output.y(neuron_idx,:)));
    nTrials             = nTotalSpikeBins/params.nTimebins;
    hist_mat            = reshape(output.y(neuron_idx,1:nTotalSpikeBins),params.nTimebins,nTrials);
    hist_model          = reshape(output.y_hat(neuron_idx,iM,1:nTotalSpikeBins),params.nTimebins,nTrials);
    
    %Show for these trial types the avg firing rate and model estimated firing rate:
    splitidxs       = {};
    
    splitidxs{1}    = strcmp(temptrialData.trialType,'P') & ismember(temptrialData.vecResponse,[1 2]);
    splitidxs{2}    = strcmp(temptrialData.trialType,'P') & temptrialData.vecResponse==3;

    for iTrial = 1:length(temptrialData.session_ID)
        shiftbins = round(temptrialData.responseLatency(iTrial) / params.binsize);
        if ~isnan(shiftbins)
            hist_mat(:,iTrial) = [hist_mat(shiftbins+1:end,iTrial); NaN(shiftbins,1)];
            hist_model(:,iTrial) = [hist_model(shiftbins+1:end,iTrial); NaN(shiftbins,1)];
        end
    end
        
    nSplits = length(splitidxs);
    
    hist_mat_mean           = NaN(nSplits,params.nTimebins);
    hist_mat_model_mean     = NaN(nSplits,params.nTimebins);
    
    for iSplit = 1:nSplits
        %Get the mean over trials::
        hist_mat_mean(iSplit,:)             = nanmean(hist_mat(:,splitidxs{iSplit}),2);
        hist_mat_model_mean(iSplit,:)       = nanmean(hist_model(:,splitidxs{iSplit}),2);
    end
    
    
    iM = 2; %show only for full model:
        
%     subplot(ceil(sqrt(nExNeurons)),floor(sqrt(nExNeurons)),iNeuron); hold all;
    subplot(nExNeurons,2,(iNeuron-1)*2+1); hold all;
%     subplot(nExNeurons,1,iNeuron); hold all;
    CurrLineHandles = [];
    for iSplit = 1:2
        
        CurrLineHandles(end+1) = plot(params.xtime,hist_mat_mean(iSplit,:),'-','LineWidth',2,'Color',colors{iSplit}); %#ok<SAGROW>
        plot(params.xtime,hist_mat_model_mean(iSplit,:),':','LineWidth',2,'Color',colors{iSplit});
        xlim([-5e5 5e5]);
        ylabel(''); xlabel('');
        set(gca,'XTick',[],'YTick',[]);
        plot([0 0],[0 100],'k:','LineWidth',2)
        
        ylim([0 ceil(max(hist_mat_mean(:)))])
        set(gca,'YTick',ceil(max(hist_mat_mean(:))))
    end
        %     Figure makeup:
%         xlabel('Time (sec)','FontSize',15)    %Label x-axis
%         ylabel('Hz','FontSize',15)            %Label y-axis
%         set(gca, 'XTick', params.xtime(1:500e3/params.binsize:params.nTimebins), 'XTickLabels', params.xtime(1:500e3/params.binsize:params.nTimebins)/1e6,'FontSize', 20)
        set(gca, 'FontSize', 20);
        legend(CurrLineHandles,labels,'FontSize',15);
        legend(gca,'boxoff');
    end
end
tightfig(gcf);

%% Make figure of explained variance;
iExp = 3;
% iExp = 1;

expidx                          = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

datatoplot                      = var_expl_splits(expidx,:);
meantoplot                      = nanmean(var_expl_splits(expidx,:),1);
errortoplot                     = nanstd(var_expl_splits(expidx,:),1) / sqrt(sum(expidx));

figure; hold all; set(gcf,'units','normalized','Position',[0.09 0.45 0.2 0.35],'color','w');

h = [];
for i = 2:params.nSplits %Give bars colors:
    h(end+1) = bar(i,meantoplot(i),0.8); %Plot bars
%     h.FaceColor = params.colorsplits{i};
    set(h(end),'FaceColor',params.colorsplits{i});
    
    errorbar(i, meantoplot(i), errortoplot(i), 'k.'); %show errorbar
end

set(gca,'XTick',2:params.nSplits,'XTickLabel',strrep(params.varsplits(2:end),'var_',' '),'XTickLabelRotation',45,'Fontsize',15);
ylabel('Explained variance')
xlim([1.5 params.nSplits+0.5])
ylim([0 0.03])
nNeurons = sum(expidx);
variance = reshape(datatoplot(:,2:end),(params.nSplits-1)*nNeurons,1);
groups = reshape(repmat(1:(params.nSplits-1),nNeurons,1),(params.nSplits-1)*nNeurons,1);
tbl = table(variance,groups);

tempBF = bf.anova(tbl,'variance~groups');

%% 
for i = 2:params.nSplits %Give bars colors:
    fprintf('%s (%2.1f%%)\n',strrep(params.varsplits{i},'var_',' '),meantoplot(i)/sum(meantoplot(2:end))*100)
end

%% Make figure of explained variance split by cell type:
iExp = 3;
% iExp = 1;

expidx                          = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

figure; hold all; set(gcf,'units','normalized','Position',[0.09 0.45 0.2 0.35],'color','w');

h = [];
for iSplit = 2:params.nSplits %Give bars colors:
    for iType = 1:2
        idx                             = expidx & spikeData.celltype==iType;

        tempmean                      = nanmean(var_expl_splits(idx,iSplit),1);
        temperror                     = nanstd(var_expl_splits(idx,iSplit),1) / sqrt(sum(idx));

        h(end+1) = bar(iSplit+(iType-1)*0.5,tempmean,0.4); %Plot bars
        tempcolor = params.colorsplits{iSplit};
        if iType ==2
            tempcolor = min([1 1 1; tempcolor + 0.6],[],1);
        end
        set(h(end),'FaceColor',tempcolor);
        errorbar(iSplit+(iType-1)*0.5, tempmean, temperror, 'k.'); %show errorbar
    end
end

set(gca,'XTick',2:params.nSplits,'XTickLabel',strrep(params.varsplits(2:end),'var_',' '),'XTickLabelRotation',45,'Fontsize',15);
ylabel('Explained variance')
xlim([1.5 params.nSplits+1])
ylim([0 0.08])

% nNeurons = sum(idx);
% variance = reshape(datatoplot(:,2:end),(params.nSplits-1)*nNeurons,1);
% groups = reshape(repmat(1:(params.nSplits-1),nNeurons,1),(params.nSplits-1)*nNeurons,1);
% tbl = table(variance,groups);
% tempBF = bf.anova(tbl,'variance~groups');

%% Make figure;
FullModelIdx    = strcmp(params.modelString,'Full');

iExp                                = 3;
expidx                              = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
expidx                              = expidx & spikeData.celltype==1;

weights_splitmat = NaN(params.nNeurons,params.nSplits);

templabels                                      = ['Offset'; output.x_label(output.x_modelidx(FullModelIdx,:))];

for iS = 1:params.nSplits
    var_idx                                     = ismember(templabels,params.(params.varsplits{iS}))';
    weights                                     = squeeze(output.weights(:,FullModelIdx,var_idx));
    weights_splitmat(:,iS)                      = nansum(weights,2);
    %         weights_splitmat(iExp,iS,1:size(weights,1)) = mean(weights,2);
end

datatoplot                      = weights_splitmat(expidx,:);
meantoplot                      = nanmean(weights_splitmat(expidx,:),1);
errortoplot                     = nanstd(weights_splitmat(expidx,:),1) / sqrt(sum(expidx));

figure; hold all; set(gcf,'units','normalized','Position',[0.09 0.45 0.2 0.35],'color','w');

h = [];
for i = 2:params.nSplits %Give bars colors:
    h(end+1) = bar(i,meantoplot(i),0.8); %Plot bars
%     h.FaceColor = params.colorsplits{i};
    set(h(end),'FaceColor',params.colorsplits{i});
    
    errorbar(i, meantoplot(i), errortoplot(i), 'k.'); %show errorbar
end

set(gca,'XTick',2:params.nSplits,'XTickLabel',strrep(params.varsplits(2:end),'var_',' '),'XTickLabelRotation',45,'Fontsize',15);
ylabel('Sum of Weights')
% xlim([1.5 5.5])

%% Explained variance over depth:
%Histogram with binning on depth:
binedges                = -100:200:1400;
binticks                = binedges(1:end-1)+50;
nBins                   = length(binedges)-1;

iExp                    = 3;
expidx                  = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
% expidx                              = true(size(spikeData.session_ID));

laminardepthfig = figure; set(gcf,'units','normalized','Position',[0.05 0.5 0.17 0.4],'color','w');
hold all; set(laminardepthfig, 'DefaultLineLineWidth', 2);

EVdepthmap         = zeros(params.nSplits,nBins); %Init variable
EVdepthmap_err         = zeros(params.nSplits,nBins); %Init variable

for iSplit = 2:params.nSplits %Give bars colors:

    for iBin = 1:nBins
        idx                         = expidx & spikeData.ChannelY>=binedges(iBin) & spikeData.ChannelY<binedges(iBin+1);
        temp                        = var_expl_splits(idx,iSplit);
        if sum(idx)>=2
            EVdepthmap(iSplit,iBin)         = nanmean(temp);
            EVdepthmap_err(iSplit,iBin)     = nanstd(temp) / sqrt(sum(idx));
        end
    end
    
    EVdepthmap(isnan(EVdepthmap)) = 0; %set to zero if no neurons sampled at this depth
    
%     EVdepthmap(iSplit,:) = EVdepthmap(iSplit,:) / max(EVdepthmap(iSplit,:));
%     plot(EVdepthmap(iSplit,:),binticks,'Color',params.colorsplits{iSplit})
    h = errorbarxy(EVdepthmap(iSplit,:),binticks,EVdepthmap_err(iSplit,:),zeros(size(binticks)),{'k.-', 'k', 'k'});
%     h = errorbarxy(EVdepthmap(iSplit,:),binticks,EVdepthmap_err(iSplit,:),zeros(size(binticks)),{'k.-', 'k', {'Color',params.colorsplits{iSplit}}});
    h.hMain.Color = params.colorsplits{iSplit};
    for iX = 1:7
        for iY = 1:6
            set(h.hErrorbar(iX,iY),'Color',params.colorsplits{iSplit});
        end
    end
    
%     'Color',params.colorsplits{iSplit})
    %Figure make up:
    ylabel('Depth (um)','FontSize', 15)
    xlabel('EV','FontSize',15)
    set(gca,'YDir','reverse','linewidth',2) %reverse depth
    Ytickies = [0 300 600 900 1200];
    set(gca,'YTick',Ytickies,'fontsize',15,'FontWeight','bold')
    set(gca,'fontsize',15,'FontWeight','bold')
    xlim([0 0.04])
    ylim([-50 1250])
%     legend boxoff
end
legend(strrep(params.varsplits(2:7),'_',' ')); legend boxoff;




%% 

% params.posthoctest  = 'bonferroni';
% params.alpha        = 0.01;
% 
% %Statistical testing:
% groups              = repmat(1:params.nExperiments*params.nSplits,500,1); groups = reshape(groups,params.nExperiments*params.nSplits*500,1);
% datatotest          = reshape(datatoplot,params.nExperiments*params.nSplits*500,1); %reshape to one column vector
% groups              = groups(~isnan(datatotest)); %filter out nans
% datatotest          = datatotest(~isnan(datatotest)); %filter out nans
% %perform kruskal wallis nonparametric anova:
% [p,table,stats]     = kruskalwallis(datatotest,groups,'off');
% comptable           = multcompare(stats,'display','off','alpha',params.alpha,'ctype',params.posthoctest); %do posthoc
% comptable           = comptable(comptable(:,end)<params.alpha,:); %Filter only significant
% sigstar(mat2cell(xpos(comptable(:,1:2)),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify

% %% Make figure;
% 
% %Init plotvar:
% datatoplot      = NaN(500,params.nExperiments,params.nSplits);
% meantoplot      = NaN(params.nExperiments,params.nSplits);
% errortoplot     = NaN(params.nExperiments,params.nSplits);
% 
% for iExp = 1:params.nExperiments %For each experiment copute the mean and error by averaging expl var over neurons
%     expidx                              = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
%     datatoplot(1:sum(expidx),iExp,:)    = var_expl_splits(expidx,:);
%     meantoplot(iExp,:)                  = nanmean(var_expl_splits(expidx,:),1);
%     errortoplot(iExp,:)                 = nanstd(var_expl_splits(expidx,:),1) / sqrt(sum(expidx));
% end
% 
% figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.05 0.9 0.55],'color','w');
% h = bar(meantoplot,1); %Plot bars
% 
% for i = 1:params.nSplits %Give bars colors:
%     h(i).FaceColor = params.colorsplits{i};
% end
% 
% % Calculating the width for each bar group
% groupwidth = min(0.8, params.nSplits/(params.nSplits + 1.5));
% xpos = [];
% for i = 1:params.nSplits
%     xpos(end+1:end+3) = (1:params.nExperiments) - groupwidth/2 + (2*i-1) * groupwidth / (2*params.nSplits);
%     errorbar(xpos(end-2:end), meantoplot(:,i), errortoplot(:,i), 'k.'); %show errorbar
% end
% set(gca,'XTick',1:3,'XTickLabel',params.ExperimentLabels,'XTickLabelRotation',0,'Fontsize',10);
% legend(h,strrep(params.varsplits,'_',' ')); legend boxoff;
% ylabel('Explained variance')
% 
% params.posthoctest  = 'bonferroni';
% params.alpha        = 0.01;
% 
% %Statistical testing:
% groups              = repmat(1:params.nExperiments*params.nSplits,500,1); groups = reshape(groups,params.nExperiments*params.nSplits*500,1);
% datatotest          = reshape(datatoplot,params.nExperiments*params.nSplits*500,1); %reshape to one column vector
% groups              = groups(~isnan(datatotest)); %filter out nans
% datatotest          = datatotest(~isnan(datatotest)); %filter out nans
% %perform kruskal wallis nonparametric anova:
% [p,table,stats]     = kruskalwallis(datatotest,groups,'off');
% comptable           = multcompare(stats,'display','off','alpha',params.alpha,'ctype',params.posthoctest); %do posthoc
% comptable           = comptable(comptable(:,end)<params.alpha,:); %Filter only significant
% sigstar(mat2cell(xpos(comptable(:,1:2)),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify

% %% Figure Weights:
% FullModelIdx    = strcmp(params.modelString,'Full');
% 
% weights_splitmat = NaN(length(params.Experiments),params.nSplits,params.nNeurons);
% for iExp = 1:length(params.Experiments)
%     expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
%     
%     templabels                                      = ['Offset'; output.x_label(output.x_modelidx(FullModelIdx,:))];
%     
%     for iS = 1:params.nSplits
%         var_idx                                     = ismember(templabels,params.(params.varsplits{iS}))';
%         weights                                     = squeeze(output.weights(expidx,FullModelIdx,var_idx));
%         weights_splitmat(iExp,iS,1:size(weights,1)) = nansum(weights,2);
% %         weights_splitmat(iExp,iS,1:size(weights,1)) = mean(weights,2);
%     end
% end
% 
% meantoplot          = nanmean(abs(weights_splitmat(:,:,:)),3);
% errortoplot         = nanstd(abs(weights_splitmat(:,:,:)),[],3);
% 
% nNeurons_Exp        = sum(~isnan(weights_splitmat(:,1,:)),3);
% for iExp = 1:params.nExperiments
%     errortoplot(iExp,:)         = errortoplot(iExp,:) / sqrt(nNeurons_Exp(iExp));
% end
% 
% figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.05 0.9 0.55],'color','w');
% h = bar(meantoplot,1);
% 
% for i = 1:params.nSplits
%     h(i).FaceColor = params.colorsplits{i};
% end
% 
% % Calculating the width for each bar group
% groupwidth = min(0.8, params.nSplits/(params.nSplits + 1.5));
% for i = 1:params.nSplits
%     x = (1:params.nExperiments) - groupwidth/2 + (2*i-1) * groupwidth / (2*params.nSplits);
%     errorbar(x, meantoplot(:,i), errortoplot(:,i), 'k.');
% end
% set(gca,'XTick',1:3,'XTickLabel',params.ExperimentLabels,'XTickLabelRotation',0,'Fontsize',10);
% legend(h,strrep(params.varsplits,'_',' ')); legend boxoff;


%% Figure correlations Weights:
FullModelIdx    = strcmp(params.modelString,'Full');

weights_splitmat = NaN(length(params.Experiments),params.nSplits,params.nNeurons);
for iExp = 1:length(params.Experiments)
    expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    templabels                                      = ['Offset'; output.x_label(output.x_modelidx(FullModelIdx,:))];
    
    for iS = 1:params.nSplits
        var_idx                                     = ismember(templabels,params.(params.varsplits{iS}))';
        weights                                     = squeeze(output.weights(expidx,FullModelIdx,var_idx));
        weights_splitmat(iExp,iS,1:size(weights,1)) = nansum(weights,2);
%                 weights_splitmat(iExp,iS,1:size(weights,1)) = nansum(abs(weights),2);
%         weights(weights==0) = NaN;
%         weights_splitmat(iExp,iS,1:size(weights,1)) = nanmean(weights,2);
        
    end
end

iExp = 3;

figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.05 0.8 0.85],'color','w');
for iX = 2:params.nSplits
    for iY = 2:params.nSplits
        subplot(params.nSplits-1,params.nSplits-1,(iY-2)*(params.nSplits-1)+iX-1)
        scatter(weights_splitmat(iExp,iX,:),weights_splitmat(iExp,iY,:),15,[0.2 0.6 0.6],'filled');
        
        Xdata       = squeeze(weights_splitmat(iExp,iX,:));
        Ydata       = squeeze(weights_splitmat(iExp,iY,:));
        idx         = ~isnan(Xdata) & ~isnan(Ydata);
        Xdata       = Xdata(idx);
        Ydata       = Ydata(idx);
        
        if iX~=iY
            tempbf      = bf.corr(Xdata,Ydata);
            bfsymb      = MOL_BFtoSymbol(tempbf);
            text(0,0.5,bfsymb,'FontSize',15)
        end
        
        %Figure makeup:
        xlim([-ceil(max(abs(Xdata))*2)/2-1 ceil(max(abs(Xdata))*2)/2+1]) %set limits, round to halfs
        ylim([-ceil(max(abs(Ydata))*2)/2-1 ceil(max(abs(Ydata))*2)/2+1])

        if iY == params.nSplits
            xlabel(strrep(params.varsplits{iX},'var_',''),'FontSize',15);
        end
        if iX == 2
            ylabel(strrep(params.varsplits{iY},'var_',''),'FontSize',15);
        end
    end
end

tightfig(gcf);

%% Figure correlations variance explained:
FullModelIdx    = strcmp(params.modelString,'Full');

figure; hold all; set(gcf,'units','normalized','Position',[0.05 0.05 0.8 0.85],'color','w');
iExp = 3;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

for iX = 2:params.nSplits
    for iY = 2:params.nSplits
        subplot(params.nSplits-1,params.nSplits-1,(iY-2)*(params.nSplits-1)+iX-1)
        scatter(var_expl_splits(expidx,iX),var_expl_splits(expidx,iY),25,[0.2 0.6 0.6],'filled');

        Xdata       = squeeze(var_expl_splits(expidx,iX));
        Ydata       = squeeze(var_expl_splits(expidx,iY));
        idx         = ~isnan(Xdata) & ~isnan(Ydata);
        Xdata       = Xdata(idx);
        Ydata       = Ydata(idx);
        
        if iX~=iY
            tempbf      = bf.corr(Xdata,Ydata);
            bfsymb      = MOL_BFtoSymbol(tempbf);
            text(0,0.05,bfsymb,'FontSize',15)
        end
        
        %Figure makeup:
        xlim([-0.1 0.1]) %set limits, round to halfs
        xlim([-0.1 0.1]) %set limits, round to halfs

%         xlim([-ceil(max(abs(Xdata))*2)/2 ceil(max(abs(Xdata))*2)/2]) %set limits, round to halfs
%         ylim([-ceil(max(abs(Ydata))*2)/2 ceil(max(abs(Ydata))*2)/2])

        xlabel(strrep(params.varsplits{iX},'var_',''),'FontSize',15);
        ylabel(strrep(params.varsplits{iY},'var_',''),'FontSize',15);
    end
end

tightfig(gcf);

%% Computing significance of EV as above a certain percent:

iExp = 3;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

varselec = [2 3 5];
% varselec = [2 3 4];

for iX = 1:length(varselec)
%     sigvarmat(iX,:) = abs(var_expl_splits(expidx,varselec(iX)))>0.01;
    sigvarmat(iX,:) = var_expl_splits(expidx,varselec(iX))>0.01;
%     sigvarmat(iX,:) = abs(var_expl_splits(expidx,varselec(iX)))>0.01;
    
end

%% Computing significance of EV versus shuffling firing rate and prediction:

iExp            = 3;
expidx          = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
varselec        = [2 3 5];

nShuffles       = 1000;
nNeurons        = size(sigvarmat,2);
alphathr        = 0.001;

idx             = params.xtime>0 & params.xtime<0.5e6;
idx             = repmat(idx,1,params.nMaxTrials);

temp_var_expl_splits = var_expl_splits(expidx,:);
evshuff         = NaN(nShuffles,1);

for iX = 1:length(varselec)
    fprintf('Computing shuffled EV across neurons for variable %d\n',iX)
    for iN = 1:nNeurons
        for iS = 1:nShuffles
            
            y_hat_shuff                = squeeze(output.y_hat_split(iN,varselec(iX),:)); %take estimate for this neuron, this variable
            y_hat_shuff                = reshape(y_hat_shuff,params.nTimebins,params.nMaxTrials); %reshape to time by trial 
            
            nTempTrials                = find(isnan(y_hat_shuff(1,:)),1)-1; %find how many trials this neuron was recorded for
            y_hat_shuff                = y_hat_shuff(:,[randperm(nTempTrials) nTempTrials+1:params.nMaxTrials]); %permute these trials, keep NaNs for higher trial numbers the same
            y_hat_shuff                = reshape(y_hat_shuff,params.nTimebins*params.nMaxTrials,1); %reshape again to vector
            
%             y_hat_shuff                = squeeze(output.y_hat_split(iN,varselec(iX),:));
%             y_hat_shuff                = reshape(y_hat_shuff,params.nTimebins,nTrials);
%             y_hat_shuff                = y_hat_shuff(:,randperm(size(y_hat_shuff, 2)));
%             y_hat_shuff                = reshape(y_hat_shuff,params.nTimebins*nTrials,1);
%             
%             y_hat_shuff         = squeeze(output.y_hat_split(iN,varselec(iX),idx));
            evshuff(iS)               = 1 - nanvar(output.y(iN,idx)' - y_hat_shuff(idx)) / nanvar(output.y(iN,idx)); %compute EV based on trial-shuffled prediction:

            %Compute variance explained by this split (subselection of predictors) for this neuron (across certain time range)
%             var_expl_splits(iNeuron,iS)     = 1 - nanvar(output.y(iNeuron,idx)' - squeeze(output.y_hat_split(iNeuron,iS,idx))) / nanvar(output.y(iNeuron,idx));
            
        end
        sigvarmat(iX,iN) = temp_var_expl_splits(iN,varselec(iX)) > prctile(evshuff,(1-alphathr)*100);

    end
    
%      = abs(var_expl_splits(expidx,varselec(iX)))>0.01;
    
end



%%
params.labels_venn      = {'Vis' 'Aud' 'Hit' 'Vis-Aud' 'Vis-Hit' 'Aud-Hit' 'All'};

frac_sign       = NaN(7,1); %store the fraction of neurons significantly coding for this variable

%compute fraction of overlap for each combination:
frac_sign(1)    = sum(sigvarmat(1,:) & ~sigvarmat(2,:) & ~sigvarmat(3,:)) / sum(expidx);
frac_sign(2)    = sum(~sigvarmat(1,:) & sigvarmat(2,:) & ~sigvarmat(3,:)) / sum(expidx);
frac_sign(3)    = sum(~sigvarmat(1,:) & ~sigvarmat(2,:) & sigvarmat(3,:)) / sum(expidx);

frac_sign(4)    = sum(sigvarmat(1,:) & sigvarmat(2,:) & ~sigvarmat(3,:)) / sum(expidx);
frac_sign(5)    = sum(sigvarmat(1,:) & ~sigvarmat(2,:) & sigvarmat(3,:)) / sum(expidx);
frac_sign(6)    = sum(~sigvarmat(1,:) & sigvarmat(2,:) & sigvarmat(3,:)) / sum(expidx);
frac_sign(7)    = sum(sigvarmat(1,:) & sigvarmat(2,:) & sigvarmat(3,:)) / sum(expidx); %triple combo

frac_noresp     = sum(~any(sigvarmat,1)) / sum(expidx);

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.3 0.3 0.3 0.4],'color','w')
%     [H,S] = venn(frac_sign,'FaceColor',{'r','y','b'},'FaceAlpha', 0.5,'EdgeColor','black');
[H,S] = venn(frac_sign,'FaceColor',params.colorsplits(varselec),'FaceAlpha',0.3,'EdgeColor',params.colorsplits(varselec));

%Now label each zone:
for i = 1:7
%     textstring = sprintf('%s %2.0f%%',params.labels_venn{i},frac_sign(i)*100);
    textstring = sprintf('%2.0f%%',frac_sign(i)*100);
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),textstring)
end
textstring = sprintf('%s %2.0f%%','Non-responsive',frac_noresp*100);
text(0.4,0.74,textstring)

%% 
params.example_cell_IDs = {
    '20102018081321110' %1 Visual cell
    '20112018081021033' %2 Auditory cell
    '20092018082321319' %3 MS cell
    '20302020011611111' %4 MS cell
    '20122018081421187' %5 Hit cell
    '20122018081421163' %6 Hit cell
    '20312020011511420' %7 Multiplex cell
    '20352019121921873' %8 lick cell
    '20092018082321322' %9 outcome history cell
    '20092018082321324' %10 choice history cell
    };

%Now label each zone:
for i = 1:10
    idx = strcmp(spikeData.cell_ID(expidx),params.example_cell_IDs{i});
    if any(idx)
        venn_cat = [];
        if any(sigvarmat(:,idx),1)
            if sigvarmat(1,idx) && ~sigvarmat(2,idx) && ~sigvarmat(3,idx); venn_cat = 1;
            elseif ~sigvarmat(1,idx) && sigvarmat(2,idx) && ~sigvarmat(3,idx); venn_cat = 2;
            elseif ~sigvarmat(1,idx) && ~sigvarmat(2,idx) && sigvarmat(3,idx); venn_cat = 3;
            elseif sigvarmat(1,idx) && sigvarmat(2,idx) && ~sigvarmat(3,idx); venn_cat = 4;
            elseif sigvarmat(1,idx) && ~sigvarmat(2,idx) && sigvarmat(3,idx); venn_cat = 5;
            elseif ~sigvarmat(1,idx) && sigvarmat(2,idx) && sigvarmat(3,idx); venn_cat = 6;
            elseif sigvarmat(1,idx) && sigvarmat(2,idx) && sigvarmat(3,idx); venn_cat = 7;
            end
            
            %     textstring = sprintf('%s %2.0f%%',params.labels_venn{i},frac_sign(i)*100);
            textstring = sprintf('%1.0f%',i);
            text(S.ZoneCentroid(venn_cat,1), S.ZoneCentroid(venn_cat,2)-0.01*i,textstring,'Color',[0.3 0.3 0.3])
        else fprintf('Example number %d not in venn diagrams\n',i)
        end
    end
end
textstring = sprintf('%s %2.0f%%','Non-responsive',frac_noresp*100);
text(0.4,0.74,textstring)

%% 

nShuffles       = 1000;
nNeurons        = size(sigvarmat,2);
tempedges       = 0:0.005:0.25;

frac_sign_shuffle       = NaN(7,nShuffles); %store the fraction of neurons significantly coding for this variable

for iS = 1:nShuffles
    for i = 1:3
        temp                = sigvarmat(i,:); 
        sigvarshuff(i,:)    = temp(randperm(nNeurons));
    end

    %compute fraction of overlap for each combination:
    frac_sign_shuffle(1,iS)    = sum(sigvarshuff(1,:) & ~sigvarshuff(2,:) & ~sigvarshuff(3,:)) / sum(expidx);
    frac_sign_shuffle(2,iS)    = sum(~sigvarshuff(1,:) & sigvarshuff(2,:) & ~sigvarshuff(3,:)) / sum(expidx);
    frac_sign_shuffle(3,iS)    = sum(~sigvarshuff(1,:) & ~sigvarshuff(2,:) & sigvarshuff(3,:)) / sum(expidx);
    
    frac_sign_shuffle(4,iS)    = sum(sigvarshuff(1,:) & sigvarshuff(2,:) & ~sigvarshuff(3,:)) / sum(expidx);
    frac_sign_shuffle(5,iS)    = sum(sigvarshuff(1,:) & ~sigvarshuff(2,:) & sigvarshuff(3,:)) / sum(expidx);
    frac_sign_shuffle(6,iS)    = sum(~sigvarshuff(1,:) & sigvarshuff(2,:) & sigvarshuff(3,:)) / sum(expidx);
    frac_sign_shuffle(7,iS)    = sum(sigvarshuff(1,:) & sigvarshuff(2,:) & sigvarshuff(3,:)) / sum(expidx); %triple combo
    
end

figure; set(gcf,'units','normalized','Position',[0.2 0.03 0.18 0.88],'color','w')
for i = 1:7
    subplot(7,1,i); hold all;
    title(params.labels_venn{i})
    plot(tempedges(1:end-1),histcounts(frac_sign_shuffle(i,:),tempedges,'normalization','probability'),'Color',[0.7 0.6 0.4]);
%     plot_shaded(tempedges(1:end-1),histcounts(frac_sign_shuffle(i,:),tempedges,'normalization','probability'),'Alpha',0,'Color',[0.7 0.6 0.4]);
    hold all;
    if frac_sign(i)>prctile(frac_sign_shuffle(i,:),97.5)
        plot([frac_sign(i) frac_sign(i)],[0 0.2],'r')
    elseif frac_sign(i)<prctile(frac_sign_shuffle(i,:),2.5)
        plot([frac_sign(i) frac_sign(i)],[0 0.2],'b')
    else
        plot([frac_sign(i) frac_sign(i)],[0 0.2],'k')
    end
    xlim([0 0.25])
    if i~=7
        set(gca,'YTick',[],'XTickLabel',[],'XTick',[0 0.05 0.1 0.15 0.2 0.25])
    else set(gca,'YTick',[],'XTick',[0 0.05 0.1 0.15 0.2 0.25])
        xlabel('Proportion')
    end
    %     frac_sign(i)>frac_sign_shuffle
    
end


%% Figure
iExp = 3;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

varselec = [2 3 5];
% varselec = [2 3 4];
% varselec = [2 6 7];
figure; set(gcf,'units','normalized','Position',[0.3 0.3 0.3 0.4],'color','w')
% scatter3(var_expl_splits(expidx,varselec(1)),var_expl_splits(expidx,varselec(2)),var_expl_splits(expidx,varselec(3)),25,[0.7 0.3 0.4],'filled');
% C = zeros(sum(expidx),3);

temp = var_expl_splits(expidx,varselec);

% temp = log10(var_expl_splits(expidx,varselec));

% C = ones(sum(expidx),3) - temp(:,1)/max(temp);
colortemp = temp;
colortemp(colortemp>0.1) = 0.1;
colortemp(colortemp<0.001) = 0.001;

C = log10(colortemp(:,[2 3 1]));
C = C - repmat(min(C),sum(expidx),1);
C = C ./ repmat(max(C),sum(expidx),1);
% C = 1 - C;
C = C.^1.4; %make everything slightly darker to avoid too yellow cells that have EV along three axes

scatter3(temp(:,1),temp(:,2),temp(:,3),35,C,'filled');

xlabel('Visual')
ylabel('Auditory')
zlabel('Hit/miss')
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
xlim([0.0001 0.1])
ylim([0.0001 0.1])
zlim([0.0001 0.4])

grid off
grid on
set(gca,'XTick',[0.001 0.01 0.1],'YTick',[0.001 0.01 0.1],'ZTick',[0.001 0.01 0.1])
set(gca,'XTickLabels',[0.001 0.01 0.1],'YTickLabels',[0.001 0.01 0.1],'ZTickLabels',[0.001 0.01 0.1])
set(gca,'XMinorTick','off','YMinorTick','off','ZMinorTick','off')

view(19.3,14)


%% 

iExp = 3;
expidx = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));

% varselec = 2:5;
varselec = [2 3 5];
% varselec = [2 6 7];
nVars = length(varselec);
for iX = 1:nVars
    sigvarmat(iX,:) = abs(var_expl_splits(expidx,varselec(iX)))>0.01;
end

frac_sign       = NaN(nVars); %store the fraction of neurons significantly coding for this variable

%Loop over variables and compute fraction of total neurons
for iX = 1:nVars
    for iY = 1:nVars
        frac_sign(iX,iY) = sum(sigvarmat(iX,:) & sigvarmat(iY,:)) / sum(expidx);
    end
end

%Make figure:
figure; set(gcf,'units','normalized','Position',[0.3 0.3 0.3 0.4],'color','w')
imagesc(frac_sign);
colormap('copper')
set(gca,'XTick',1:nVars,'XTickLabel',strrep(params.varsplits(varselec),'var_',''))
set(gca,'YTick',1:nVars,'YTickLabel',strrep(params.varsplits(varselec),'var_',''),'xaxisLocation','top','XTickLabelRotation',45)

% scatter(var_expl_splits(expidx,iX),var_expl_splits(expidx,iY),25,[0.2 0.6 0.6],'filled');




