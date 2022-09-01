%%
startover


%% Supplementary figure: 
% [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},{'2009' '2010' '2011' '2012' '2013' '2019' '2020' '2021' '2022' '2023'},[],{'sessionData' 'trialData'});
[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% General settings:
params.sumTrialsCondition       = 30; %Combined minimum number of trials per opto condition to fit the model
params.minPerf                  = 0.3;
params.minPhotoStimPower        = 2;

params                          = MOL_getColors_CHDET(params);

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

sesids              = sessionData.session_ID(visperf>0.3 & auperf>0.3);
fprintf('Removed %d/%d sessions with low behavioral accuracy\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Filter sessions with optogenetic manipulation in PPC:
sesids            = unique(sessionData.session_ID(sessionData.UseOpto==1 & ...
    strcmp(sessionData.PhotostimArea,'PPC') & ...
    sessionData.Photostimpower >= params.minPhotoStimPower));
fprintf('Selected %d/%d sessions with optical manipulation in V1 or PPC\n',length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Filter out session with very poor performance 
nSessions = length(sessionData.session_ID);
for iSes = 1:nSessions
    idx_trial = strcmp(trialData.session_ID,sessionData.session_ID(iSes));
    idx_flag(iSes) = numel(unique(trialData.visualOriChangeNorm(idx_trial)))>=4 && numel(unique(trialData.audioFreqChangeNorm(idx_trial)))>=4;
end
idx_flag            = idx_flag & ~strcmp(sessionData.session_ID,'20122018080318')'; %Excl criteria, see Methods manuscript (FArate)

fprintf('Selected %d/%d sessions with sufficient trial conditions for psychometrics\n',sum(idx_flag),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sessionData.session_ID(idx_flag),sessionData,trialData);

%% Trim post change delayed inhibition (for V1 - 2nd bump project):
trialData.PostChangeOptoStart(trialData.PostChangeOptoStart<0)=0;

%% Filter out sessions with only delayed or prechange photostimulation:
sesids              = unique(trialData.session_ID(trialData.PostChangeOptoStart==0));
fprintf('Removed %d/%d sessions with prechange or delayed opto inhibition\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Initialize structure for saving output fit parameters:
nSessions               = length(sessionData.session_ID);
FullParam_ctrl          = NaN(nSessions,8); %Init output variables
FullParam_opto          = NaN(nSessions,8); %Init output variables
TotalResp               = NaN(5,5,3,3,nSessions); %Init matrix for storing all hitrate data

%% Loop over sessions:
for iSes = 1:nSessions
    fprintf('Fitting session %d/%d\n',iSes,nSessions);
    sesid                                                       = sessionData.session_ID(iSes);
    [tempsessionData,temptrialData]                             = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    
    trialFields = fieldnames(temptrialData);
    for iF = 1:length(trialFields)
        if any(strfind(trialFields{iF},'audio')) && iscell(temptrialData.(trialFields{iF}))
            temptrialData.(trialFields{iF}) = cell2vec(temptrialData.(trialFields{iF}))';
        end
    end
    
    [visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMatOpto(tempsessionData,temptrialData);
    %output is Au X Vis X response matrix (dimension 3 (response): layer 1 is fraction auditory, 2 visual, 3 no response)

    %Rewrite fullresponse table to contingency table for model fitting (without conflict trials)
    FullRespMat_ctrl        = squeeze(FullRespMat(:,:,1,:));
    ctable_ctrl             = NaN(3,3,numel(visconditions));
    ctable_ctrl(1,:,:)      = permute(FullRespMat_ctrl(2:end,1,:),[3 1 2]);              %Auditory
    ctable_ctrl(2,:,:)      = permute(FullRespMat_ctrl(1,2:end,:),[1 3 2]);              %Visual
    ctable_ctrl(3,:,:)      = repmat(permute(FullRespMat_ctrl(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
    
    FullRespMat_opto        = squeeze(FullRespMat(:,:,2,:));
    ctable_opto             = NaN(3,3,numel(visconditions));
    ctable_opto(1,:,:)      = permute(FullRespMat_opto(2:end,1,:),[3 1 2]);              %Auditory
    ctable_opto(2,:,:)      = permute(FullRespMat_opto(1,2:end,:),[1 3 2]);              %Visual
    ctable_opto(3,:,:)      = repmat(permute(FullRespMat_opto(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
    
    %Store values:
    TotalResp(:,:,:,:,iSes)     = FullRespMat;
    
    %Align visual and auditory intensities
    %visual and auditory conditions should be normalized such that
    %value of 1 corresponds to expected asymptotic dmax
    normconditions = [auconditions / max(auconditions); visconditions / max(visconditions)];
    
    %Fit psychometric 2ADC model for control trials:
    [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable_ctrl,normconditions,[]);
    %Store parameters:
    %theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    FullParam_ctrl(iSes,:)            = theta_est; %store parameters
    
    %Fit psychometric 2ADC model for opto trials:
    [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable_opto,normconditions,[]);
    %Store parameters:
    %theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    FullParam_opto(iSes,:)            = theta_est; %store parameters
end
fprintf('\n\n Finished fitting behavioral model.\n\n\n')

%% FIGURES:
uMice                 = unique(sessionData.mousename);
nMice                 = length(uMice);
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;

for iMou = 1:nMice
    
    idx_ses                     = strcmp(sessionData.mousename,uMice{iMou});
    [tempsessionData,~]         = MOL_getTempPerSes(sessionData.session_ID(idx_ses),sessionData,trialData);
    meantheta_est_ctrl          = mean(FullParam_ctrl(idx_ses,:),1);
    meantheta_est_opto          = mean(FullParam_opto(idx_ses,:),1);
    
    % Settings:
    switch tempsessionData.auChangeUnit{1}
        case 'Hz'
            params.auprobepos = 0.5;
            params.auticks = [100 4000];
            params.auticklabels = ['Probe' num2cell(params.auticks)];
            params.auxaxislabel  = 'Delta frequency (Hz)';
            params.auystatslabel = 'Auditory threshold (Hz)';
        case 'Oct'
            params.auprobepos = 0.001;
            params.auticks = [1/256 1/64 1/8 1/2];
            params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
            params.auxaxislabel  = 'Delta frequency (Oct)';
            params.auystatslabel = 'Auditory threshold (partial octave)';
    end
    
    params.visprobepos     = 0.5;
    params.visticks        = [5 30 90];
    params.vistickslabels  = ['Probe' num2cell(params.visticks)];
    params.visxaxislabel   = 'Delta orientation (Degrees)';
    params.visystatslabel  = 'Visual threshold (Degrees)';
    
    params.yticks          = [0 0.25 0.5 0.75 1];
    
    % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est_ctrl,params);
    %Audio:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',3);
    %Visual:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',3);
    
     % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est_opto,params);
    %Audio:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),':','Color',[1 0 0],'LineWidth',3);
    %Visual:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),':','Color',[0 0 1],'LineWidth',3);
end

MOL_Psy2Sided_FigMakeup(params)

%% FIGURES:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;

for iSes = 1:nSessions
%     idx_ses                     = strcmp(sessionData.mousename,uMice{iMou});
    [tempsessionData,~]         = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);
%     meantheta_est_ctrl          = mean(FullParam_ctrl(idx_ses,:),1);
%     meantheta_est_opto          = mean(FullParam_opto(idx_ses,:),1);
    
    % Settings:
    switch tempsessionData.auChangeUnit{1}
        case 'Hz'
            params.auprobepos = 0.5;
            params.auticks = [100 4000];
            params.auticklabels = ['Probe' num2cell(params.auticks)];
            params.auxaxislabel  = 'Delta frequency (Hz)';
            params.auystatslabel = 'Auditory threshold (Hz)';
        case 'Oct'
            params.auprobepos = 0.001;
            params.auticks = [1/256 1/64 1/8 1/2];
            params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
            params.auxaxislabel  = 'Delta frequency (Oct)';
            params.auystatslabel = 'Auditory threshold (partial octave)';
    end
    
    params.visprobepos     = 0.5;
    params.visticks        = [5 30 90];
    params.vistickslabels  = ['Probe' num2cell(params.visticks)];
    params.visxaxislabel   = 'Delta orientation (Degrees)';
    params.visystatslabel  = 'Visual threshold (Degrees)';
    
    params.yticks          = [0 0.25 0.5 0.75 1];
    
    % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(FullParam_ctrl(iSes,:),params);
    %Audio:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0.5 0.5],'LineWidth',2);
%     plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'-','Color',[0.5 0.5 1],'LineWidth',1);
    
    %Visual:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0.5 0.5 1],'LineWidth',2);
%     plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),'-','Color',[1 0.5 0.5],'LineWidth',1);
     % Generate contingency table from fitted parameters:
    [xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(FullParam_opto(iSes,:),params);
    %Audio:
    subplot(1,2,1); hold all;
    plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),':','Color',[1 0.5 0.5],'LineWidth',2);
%     plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',[0.5 0.5 1],'LineWidth',1);
    %Visual:
    subplot(1,2,2); hold all;
    plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),':','Color',[0.5 0.5 1],'LineWidth',2);
%     plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),':','Color',[1 0.5 0.5],'LineWidth',1);
end

meantheta_est_ctrl          = median(FullParam_ctrl,1);
meantheta_est_opto          = median(FullParam_opto,1);

% Generate contingency table from fitted parameters:
[xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est_ctrl,params);
%Audio:
subplot(1,2,1); hold all;
plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',5);
% plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),'-','Color',[1 0 0],'LineWidth',5);
%Visual:
subplot(1,2,2); hold all;
plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',5);

% Generate contingency table from fitted parameters:
[xvals_fit_au,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est_opto,params);
%Audio:
subplot(1,2,1); hold all;
plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),':','Color',[1 0 0],'LineWidth',5);
%Visual:
subplot(1,2,2); hold all;
plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),':','Color',[0 0 1],'LineWidth',5);

MOL_Psy2Sided_FigMakeup(params)

%% Show average response rate figure for photoinhibition @PPC in MST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'PPC') & ...
    ismember(sessionData.Experiment,{'BehaviorConflict' 'ChangeDetectionConflict'});

params.colors_audio_opto       = {[212 0 41] [255 122 94] [10 10 10]}; 
params.colors_audio_opto       = cellfun(@(x) x/256,params.colors_audio_opto,'UniformOutput',false);

params.colors_visual_opto       = {[40 0 150] [0 173 240] [10 10 10]}; 
params.colors_visual_opto       = cellfun(@(x) x/256,params.colors_visual_opto,'UniformOutput',false);

params.lines_audio_opto       = {'-o' '--o' ':o'}; 
params.lines_visual_opto       = {'-o' '--o' ':o'}; 

fprintf( '\nHit rates:\n')
MOL_plotOptoBehaviorPsy_Rates_PPC(params,TotalResp(:,:,:,:,idx_ses))
% fprintf( '\nDprime:\n')
% MOL_plotOptoBehavior_Dprime_PPC(params,dVis(idx_ses,:,:),dAud(idx_ses,:,:))


%% Make figures of dmax and thresholds on control or on opto trials:
fitvars         = [1 2 3 4 5 6 7 8];
% colors_splits   = {'r' 'b' 'r' 'b'};
colors_splits   = {'r' 'b' 'r' 'b' 'r' 'b' 'r' 'b'};
labels_splits   = {'D''Vis' 'D''Aud' 'nVis' 'nAud' 'Au Thr' 'Vis Thr' 'cAud' 'cVis'};

% fitvars         = [1 2 7 8];
% colors_splits   = {'r' 'b' 'r' 'b'};
% labels_splits   = {'D''' 'D''' 'Au Thr' 'Vis Thr'};

figure; set(gcf,'color','w','units','normalized','Position', [0.05 0.15 .45 .6]); hold all;

for iP = 1:length(fitvars)
    subplot(2,4,iP)
    for iSes = 1:nSessions
        plot([1 2],[FullParam_ctrl(iSes,fitvars(iP)) FullParam_opto(iSes,fitvars(iP))],'.-','Color','k','MarkerSize',30,'MarkerEdgeColor',colors_splits{iP},'LineWidth',1); hold all;
    end
    
    %Classical statistics:
%     p = signrank(FullParam_ctrl(:,fitvars(iP)),FullParam_opto(:,fitvars(iP)));
%     sigstar([1 2],p) %use sigstar function to identify signficantly different conditions
%     text(1,0.8,sprintf('p=%1.3f',p),'FontSize',15)
    
    %Bayesian statistics:
    tempbf      = bf.ttest(FullParam_ctrl(:,fitvars(iP)),FullParam_opto(:,fitvars(iP)));
    bfsymb      = MOL_BFtoSymbol(tempbf);
    text(mean([1 2]),nanmean([FullParam_ctrl(:,fitvars(iP)); FullParam_opto(:,fitvars(iP))]),bfsymb,'FontSize',15)
    tempd       = computeCohen_d(FullParam_ctrl(:,fitvars(iP)),FullParam_opto(:,fitvars(iP)),'paired');
    fprintf('Bayesian ttest (%s): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',labels_splits{iP},size(FullParam_opto(:,fitvars(iP)),1),tempd,tempbf)
%     bayesstar([1 2],tempbf);
    
    %Make up:
    set(gca, 'XTick', 1:4,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
    xlim([0.7 2.3])
    if iP <=4
        ylim([0 4])
    end
    if iP == 5
        set(gca,'Yscale','log','YTick',[10 100 500])
        ylim([0.5 500])
    end
    
    ylabel(labels_splits{iP})
end

%%
