%% MOL_Zmat_PPC_TrainedNaive
% This script plots the average auditory and average visual response of all PPC neurons in trained and in naive animals

%% Parameter settings:
params                  = params_histresponse(); %Parameters for PSTH (All time is in microseconds)

params.Experiments      = {'ChangeDetectionConflictDecor' 'ChangeDetectionConflict' }; %Which versions of the task to load data from
params.ExperimentLabels = {'NE' 'MST'}; %Labels for the different experiments

params.AlignOn          = 'stimChange';      %On which timestamp to align as t=0

params.minTrialCond     = 5;
params.cscale           = [-0.3 2.2];

params.area             = 'PPC';

params.twin_resp_start  = 0e6;
params.twin_resp_stop   = 1e6;
params.conv_sigma       = 50e3;
params.conv_win          = 'gaussian';

params                  = MOL_getColors_CHDET(params);

params.savedir          = 'E:\Documents\PhD\Figures\Project CHDET\Manuscript  - PPC\Figure 1 - Task Active vs Passive PPC';

%% Get data:
[Data] = MOL_GetData('E:','CHDET',params.Experiments,{},[],{'sessionData' 'trialData_newtrials' 'spikeData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;
spikeData       = Data.spikeData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter neurons on area:
spikeFields     = fieldnames(spikeData);
idx             = ismember(spikeData.area,params.area);
for iF = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(idx);
end
fprintf('Filtered neurons based on area\n');

%% Filter out neurons based on quality:
spikeData       = MOL_filterNeurons(sessionData,trialData,spikeData);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions  with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out trials with photostimulation:
sesids              = sessionData.session_ID(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC')));
idx                 = trialData.hasphotostim==1 & ismember(trialData.session_ID,sesids);
trialFields         = fieldnames(trialData);
fprintf('Removed %d/%d trials  with photostim in PPC or V1\n',sum(idx),length(trialData.session_ID));
for iF = 1:length(trialFields)
    trialData.(trialFields{iF}) = trialData.(trialFields{iF})(~idx);
end

%% Filter out latest recording animals with passive sessions and muscimol:
sesids              = sessionData.session_ID(~ismember(sessionData.mousename,{'2044' '2045'}));
fprintf('Removed %d/%d sessions from animals 2044 and 2045 with passive epochs \n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Filter out multisensory sessions that have too low visual or auditory performance:
nSessions           = length(sessionData.session_ID);
visperf             = NaN(nSessions,1);
auperf              = NaN(nSessions,1);
for iSes = 1:nSessions
    sesidx          = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    vistrialidx     = strcmp(trialData.trialType,'X') & trialData.visualOriChangeNorm==3 & trialData.hasphotostim~=1 & sesidx;
    visperf(iSes)   = sum(trialData.correctResponse(vistrialidx)) / sum(vistrialidx);
    autrialidx      = strcmp(trialData.trialType,'Y') & trialData.audioFreqChangeNorm==3 & trialData.hasphotostim~=1 & sesidx;
    auperf(iSes)    = sum(trialData.correctResponse(autrialidx)) / sum(autrialidx);
end

sesids              = sessionData.session_ID(~(strcmp(sessionData.Experiment,'ChangeDetectionConflict') & (visperf<0.3 | auperf<0.3)));
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData, spikeData]        = MOL_getTempPerSes(sesids,sessionData,trialData,spikeData);

%% Get only sessions with neurons recorded:
[sessionData,trialData,spikeData] = MOL_getTempPerSes(unique(spikeData.session_ID),sessionData,trialData,spikeData);
fprintf('Dataset: %d mice, %d sessions, %d trials, %d neurons\n',length(unique(sessionData.mousename)),length(sessionData.session_ID),length(trialData.session_ID),length(spikeData.session_ID));

%% Main loop to get psth matrix:
params.nSplits          = 2;

params.labels_splits    = {'Visual' 'Auditory'};
params.colors_splits    = params.colors_trialtypes([2 1]);

nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing average Z-scored response for neuron        \n');
snakemat                = NaN(nNeurons,params.nTimebins,params.nSplits);

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+1));
    fprintf('%d/%d',iNeuron,nNeurons);
    
    if ~strcmp(lastsesid,spikeData.session_ID(iNeuron)) %construct new predictor matrix if neuron comes from a new session:
        %Get the relevant data for each session individually:
        temptrialData        = MOL_getTempPerSes(spikeData.session_ID(iNeuron),trialData);
        lastsesid            = spikeData.session_ID(iNeuron); %save this session_ID
    end
    
    %Compute histogram:
    events_ts               = temptrialData.(params.AlignOn);
    hist_mat                = calc_psth(events_ts,spikeData.ts{iNeuron},params);    %Construct histogram matrix:
    
    splits                  = {};
    splits{1}               = ismember(temptrialData.trialType,{'X'}); %Visual trials
    splits{2}               = ismember(temptrialData.trialType,{'Y'}); %Auditory trials
    
    for iSplit = 1:params.nSplits
        if sum(splits{iSplit})>=params.minTrialCond
            snakemat(iNeuron,:,iSplit) = mean(hist_mat(splits{iSplit},:),1);
        end
    end
end
fprintf('\n');

%% Sort snakemat
params.twin_resp_start  = 0e6;
params.twin_resp_stop   = 0.5e6;

meanresp                    = mean(snakemat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,1),2) + mean(snakemat(:,params.xtime>params.twin_resp_start & params.xtime<=params.twin_resp_stop,2),2);

[~,sortidx]                 = sort(meanresp,1,'descend');

snakemat        = snakemat(sortidx,:,:);
spikeFields     = fieldnames(spikeData);
for iF          = 1:length(spikeFields)
    spikeData.(spikeFields{iF}) = spikeData.(spikeFields{iF})(sortidx,:);
end

%% Make figure:
% params.colormap = 'parula';

params.cscale = [-0 1.5];
params.cscale = [-1.2 1.2];

params.colors_splits = {[0.2 0.2 1] [1 0.2 0.2]};

nExperiments = length(params.Experiments);
figure; set(gcf,'units','normalized','Position',[0.3 0.5 0.18 0.27],'color','w')

for iExp = 1:nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
    for iSplit = 1:length(splits)
        ax = subplot(2,2,(iSplit-1)*2+iExp);
        %  ax = subplot(2,2,(iExp-1)*2+iSplit);
        title(params.ExperimentLabels{iExp})
        imagesc(snakemat(idx_exp,:,iSplit),params.cscale); hold on;
        %         imagesc(snakemat(idx_exp & ~idx_clip,:,iSplit),params.cscale); hold on;
        plot([find(params.xtime == 0) find(params.xtime == 0)], [0 size(snakemat(idx_exp,:,iSplit),1)+0.5],':','Color',[0.6 0.6 0.6],'LineWidth',2);
        
        set(gca, 'XTick', 1:1000:length(params.xtime), 'XTickLabels', params.xtime(1:1000:length(params.xtime))/1e6,'FontSize', 20)
        set(gca, 'YTick', [1 sum(idx_exp)], 'YTickLabels', [1 sum(idx_exp)],'FontSize', 7)
        xlim([find(params.xtime == -0.8e6) find(params.xtime == 1.4e6)]);
        xlabel('Time (s)','FontSize', 7)
        if iExp==1
            ylabel('Neuron','FontSize', 7)
        end
        
        
        params.colormap = [linspace(0,params.colors_splits{iSplit}(1))' linspace(0,params.colors_splits{iSplit}(2))' linspace(0,params.colors_splits{iSplit}(3))'];
%         params.colormap = [linspace(1,params.colors_splits{iSplit}(1))' linspace(1,params.colors_splits{iSplit}(2))' linspace(1,params.colors_splits{iSplit}(3))'];
        
        temp1 = [linspace(0.3,0)' linspace(1,0)' linspace(1,0)'];
        temp2 = [linspace(0,params.colors_splits{iSplit}(1))' linspace(0,params.colors_splits{iSplit}(2))' linspace(0,params.colors_splits{iSplit}(3))'];
        params.colormap = [temp1; temp2];
        
        temp1 = [linspace(0,1)' linspace(0,1)' linspace(0,1)'];
        temp2 = ones(20,3);
        temp3 = [linspace(1,params.colors_splits{iSplit}(1))' linspace(1,params.colors_splits{iSplit}(2))' linspace(1,params.colors_splits{iSplit}(3))'];
        params.colormap = [temp1; temp2; temp3;];
        
%         params.colormap = getPyPlot_cMap('seismic',1000);
        params.colormap = getPyPlot_cMap('PuOr',1000);

        colormap(ax,params.colormap);
        if iExp==2
            c = colorbar;
            c.Position = c.Position + 1e-10;
            c.Label.String = 'Z-scored firing rate';
        end
    end
end

%% Make figure of the mean:
nExperiments = length(params.Experiments);
for iExp = 1:nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    figure; set(gcf,'units','normalized','Position',[0.4 0.3 0.32 0.38],'color','w')
    
    handles = [];
    for iSplit = 1:length(splits)
        meantoplot = nanmean(snakemat(idx_exp,:,iSplit),1);
        %         meantoplot = nanmean(abs(snakemat(idx_exp,:,iSplit)),1);
        errortoplot = nanstd(snakemat(idx_exp,:,iSplit),1)/sqrt(size(snakemat(idx_exp,:,iSplit),1));
        h = shadedErrorBar(params.xtime,meantoplot,errortoplot,{'-k','markerfacecolor',params.colors_splits{iSplit},'LineWidth',3},1);
        h.mainLine.Color = params.colors_splits{iSplit};    h.patch.FaceColor = params.colors_splits{iSplit};
        delete(h.edge(1)); delete(h.edge(2));
        handles(iSplit) = h.mainLine; hold all; %#ok<SAGROW>
    end
    set(gca, 'XTick', params.xtime(1:500:length(params.xtime)), 'XTickLabels', params.xtime(1:500:length(params.xtime))/1e6,'FontSize', 20)
    xlim([-0.3e6 1.5e6]);
    ylim([-0.2 0.7])
    %     ylim([-0.2 0.9])
    title(params.ExperimentLabels{iExp})
    ylabel('Z-scored firing rate','FontSize', 20)
    xlabel('Time (s)','FontSize', 20)
    legend(handles,params.labels_splits); legend boxoff
end

%% Compute significantly responsive units:
[sessionData,trialData,spikeData] = MOL_calc_sign_AV_PPC(sessionData,trialData,spikeData);

%% Make pie chart of the number of significantly responding neurons in each experiment:

labels_resp = {'NonResp' 'Decr' 'Incr'};
labels_mod = {'Auditory' 'Visual'};
fracmat = NaN(2,2,3);

params.colors_mods{1} = params.colormap(100,:);
params.colors_mods{2} = params.colormap(end-100,:);

nExperiments = length(params.Experiments);
figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.45 0.6],'color','w')
for iExp = 1:nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    
%     if iExp ==2
%         idx_side    = ismember(spikeData.session_ID,sessionData.session_ID(sessionData.VisualLeftCorrectSide==1));
%         idx_exp     = idx_exp & idx_side;
%     end
    
    %     figure; set(gcf,'units','normalized','Positi0n',[0.2 0.3+(iExp-1)*0.3 0.2 0.3],'color','w')
    subplot(2,2,iExp);
%     subplot(2,2,(iExp-1)*2+1)
    
    X = [sum(~(spikeData.sign_incr_aud(idx_exp) | spikeData.sign_decr_aud(idx_exp))) sum(spikeData.sign_decr_aud(idx_exp)) sum(spikeData.sign_incr_aud(idx_exp))];
    fprintf('%s: %s: %s:\n',params.ExperimentLabels{iExp},labels_mod{1},labels_resp{:})
    fprintf('%2.3f%%\n',X/sum(idx_exp))
    
    p = pie(X,labels_resp);
    p(1).FaceColor = [0.8 0.8 0.8];
    p(3).FaceColor = params.colors_mods{1};
    p(5).FaceColor = params.colors_mods{2};
    title(sprintf('%s - %s',labels_mod{1},params.ExperimentLabels{iExp}))
    fracmat(1,iExp,:) = X;
    %     saveas(gcf,fullfile(params.savedir,sprintf('Pie_PPC_AudResponsive_%s.bmp',params.ExperimentLabels{iExp})));
    %     saveas(gcf,fullfile(params.savedir,sprintf('Pie_PPC_AudResponsive_%s.eps',params.ExperimentLabels{iExp})));
    
    %     figure; set(gcf,'units','normalized','Position',[0.4 0.3+(iExp-1)*0.3 0.2 0.3],'color','w')
%     subplot(2,2,(iExp-1)*2+2)
    subplot(2,2,2+iExp);
    X = [sum(~(spikeData.sign_incr_vis(idx_exp) | spikeData.sign_decr_vis(idx_exp))) sum(spikeData.sign_decr_vis(idx_exp)) sum(spikeData.sign_incr_vis(idx_exp))];
    fprintf('%s: %s: %s:\n',params.ExperimentLabels{iExp},labels_mod{2},labels_resp{:})
    fprintf('%2.3f%%\n',X/sum(idx_exp))
    p = pie(X,labels_resp);
    p(1).FaceColor = [0.8 0.8 0.8];
    p(3).FaceColor = params.colormap(100,:);
    p(5).FaceColor = params.colormap(end-100,:);
    
%     p(1).FaceColor = [0.5 0.5 0.5];
%     p(3).FaceColor = [0 0 0.8];
%     p(5).FaceColor = [0.2 0.2 1];
    %     saveas(gcf,fullfile(params.savedir,sprintf('Pie_PPC_VisResponsive_%s.bmp',params.ExperimentLabels{iExp})));
    %     saveas(gcf,fullfile(params.savedir,sprintf('Pie_PPC_VisResponsive_%s.eps',params.ExperimentLabels{iExp})));
    fracmat(2,iExp,:) = X;
    title(sprintf('%s - %s',labels_mod{2},params.ExperimentLabels{iExp}))

    
end

for iMod = 1:2
    for iDir = 1:2
        
        n1 = fracmat(iMod,1,iDir+1); N1 = sum(fracmat(iMod,1,:));
        n2 = fracmat(iMod,2,iDir+1); N2 = sum(fracmat(iMod,2,:));
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [ones(n1,1); repmat(2,N1-n1,1); ones(n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,pval] = crosstab(x1,x2);
        
        fprintf('Chi square two samples test for %s %s, %1.6f\n',labels_mod{iMod},labels_resp{iDir+1},pval)
        
    end
end

%% Make figure of the mean:
params.idx_splits       = [spikeData.sign_incr_vis spikeData.sign_decr_vis spikeData.sign_incr_aud spikeData.sign_decr_aud];
params.idx_splits       = cat(3,[spikeData.sign_decr_vis spikeData.sign_incr_vis],[spikeData.sign_decr_aud spikeData.sign_incr_aud]);

% idx_time = params.xtime>0 & params.xtime<0.5e6;
% params.idx_splits       = [nanmean(snakemat(:,idx_time,iMod),2)>0 nanmean(snakemat(:,idx_time,iMod),2)<0 nanmean(snakemat(:,idx_time,iMod),2)>0 nanmean(snakemat(:,idx_time,iMod),2)<0];

params.lines_exps       = {':' '-'};

nExperiments = length(params.Experiments);
figure; set(gcf,'units','normalized','Position',[0.15 0.3 0.5 0.35],'color','w')
plot_time_idx = params.xtime>-0.2e6 & params.xtime<1e6;

for iMod = 1:2
    subplot(1,2,iMod); hold all;
    handles = [];
    
    for iExp = 1:nExperiments
        idx_exp     = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
%         if iExp ==2
%             idx_side    = ismember(spikeData.session_ID,sessionData.session_ID(sessionData.VisualLeftCorrectSide==1));
%             idx_exp     = idx_exp & idx_side;
%         end
        
        for iSplit = 1:2
%             idx             = idx_exp & params.idx_splits(:,iSplit+iMod-1);
            idx             = idx_exp & params.idx_splits(:,iSplit,iMod);
            
            idx_exp_tr = ismember(trialData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
            if iMod==1
                resplat                 = nanmedian(trialData.responseLatency(ismember(trialData.trialType,{'X'}) & idx_exp_tr & trialData.vecResponse==2)); %Visual trials
            elseif iMod ==2
                resplat                 = nanmedian(trialData.responseLatency(ismember(trialData.trialType,{'Y'}) & idx_exp_tr & trialData.vecResponse==1)); %Visual trials
            end
            
            meantoplot      = nanmean(snakemat(idx,:,iMod),1);
            errortoplot     = nanstd(snakemat(idx,:,iMod),1)/sqrt(size(snakemat(idx,:,iMod),1));
            
            h = shadedErrorBar(params.xtime(plot_time_idx),meantoplot(plot_time_idx),errortoplot(plot_time_idx),{params.lines_exps{iExp},'markerfacecolor',params.colors_mods{iSplit},'LineWidth',3},0);
            h.mainLine.Color = params.colors_mods{iSplit};    h.patch.FaceColor = [0.6 0.6 0.6];
            delete(h.edge(1)); delete(h.edge(2));
            if iExp==2 && iSplit==2
                plot(resplat,1.5,'v','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerSize',7)
            end

        end
        handles(iExp) = h.mainLine; hold all; %#ok<SAGROW>
    end
    
    title(params.labels_splits{iMod})
    set(gca, 'XTick',[-0.2e6 0 0.5e6 1e6], 'XTickLabels', [-0.2e6 0 0.5e6 1e6]/1e6,'FontSize', 20)
    xlim([-0.2e6 1e6]);
    ylim([-0.6 1.5])
    
    if iMod == 1
        ylabel('Z-scored firing rate','FontSize', 20)
    end
    xlabel('Time (s)','FontSize', 20)
    legend(handles,params.ExperimentLabels); legend boxoff
end
MOL_prepfigAI

%% Make bar plot of contralateral versus ipsilateral number of significantly responding neurons:

labels_resp     = {'Decr' 'Incr'};
labels_mod      = {'Auditory' 'Visual'};
fracmat         = NaN(2,2,2,2); %dims: 1) experiments (2), 2) modalities (2), 3) modality-side pairing 4) increasing/decreasing (2)

params.colors_mods{1} = params.colormap(100,:);
params.colors_mods{2} = params.colormap(end-100,:);

nExperiments = length(params.Experiments);
for iExp = 1:nExperiments
    idx_exp = ismember(spikeData.session_ID,sessionData.session_ID(strcmp(sessionData.Experiment,params.Experiments(iExp))));
    for iP = 1:2
        if iExp ==2
            idx_side    = ismember(spikeData.session_ID,sessionData.session_ID(sessionData.VisualLeftCorrectSide==iP-1));
            idx_all     = idx_exp & idx_side;
        else
            idx_all     = idx_exp;
        end
        
        if iP==1
            X = [sum(spikeData.sign_decr_aud(idx_all)) sum(spikeData.sign_incr_aud(idx_all))];
            fracmat(iExp,1,1,:) = X / sum(idx_all);
            
            X = [sum(spikeData.sign_decr_vis(idx_all)) sum(spikeData.sign_incr_vis(idx_all))];
            fracmat(iExp,2,2,:) = X / sum(idx_all);
        elseif iP==2
            X = [sum(spikeData.sign_decr_aud(idx_all)) sum(spikeData.sign_incr_aud(idx_all))];
            fracmat(iExp,1,2,:) = X / sum(idx_all);
            
            X = [sum(spikeData.sign_decr_vis(idx_all)) sum(spikeData.sign_incr_vis(idx_all))];
            fracmat(iExp,2,1,:) = X / sum(idx_all);
        end
        
    end
end

figure; set(gcf,'units','normalized','Position',[0.2 0.2 0.11 0.3],'color','w'); hold all;
iP = 1;
for iMod = 1:2
    plot([1 2],squeeze(fracmat(:,iMod,iP,1)),':.','LineWidth',1,'MarkerSize',25,'Color',params.colors_modalities{iMod})
    plot([1 2],squeeze(fracmat(:,iMod,iP,2)),'.-','LineWidth',1,'MarkerSize',25,'Color',params.colors_modalities{iMod})
end

iP = 2;
for iMod = 1:2
    plot([3 4],squeeze(fracmat(:,iMod,iP,1)),':.','LineWidth',1,'MarkerSize',25,'Color',params.colors_modalities{iMod})
    plot([3 4],squeeze(fracmat(:,iMod,iP,2)),'.-','LineWidth',1,'MarkerSize',25,'Color',params.colors_modalities{iMod})
end

xlim([0.5 4.5]); ylim([0 0.4]);
set(gca,'XTick',1:4,'XTickLabels',[params.ExperimentLabels params.ExperimentLabels],'YTick',0:0.1:0.5)
ylabel('% responsive')

contramat       = fracmat(:,:,1,:);
ipsimat         = fracmat(:,:,2,:);

[tempbf] = bf.ttest(contramat(:),ipsimat(:));
fprintf('Contra versus ipsilateral responses: n=8 fractions, BF=%3.3f\n',tempbf)

bayesstar([1.5 3.5],tempbf)
legend(labels_resp); legend boxoff
