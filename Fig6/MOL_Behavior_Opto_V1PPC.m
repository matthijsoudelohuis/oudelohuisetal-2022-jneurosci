%% MOL_Behavior_Opto_V1PPC

%% Get Data:
animals = {'2003'    '2004'    '2007'    '2009'    '2010'    '2011'   '2012'    '2013'    '2019'    '2020'    '2021'    '2022'    '2023'    '2026'    '2027'    '2030'    '2031'};

[Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},animals,[],{'sessionData' 'trialData'});
% [Data] = MOL_GetData('E:','CHDET',{'ChangeDetectionConflict' 'BehaviorConflict'},{'2009' '2010' '2011' '2012' '2013' '2019' '2020' '2021' '2022' '2023'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% General settings:
params.sumTrialsCondition       = 30; %Combined minimum number of trials per opto condition to fit the model
params.minPerf                  = 0.3;
params.minPhotoStimPower        = 2;

params                          = MOL_getColors_CHDET(params);

params.posthoctest              = 'bonferroni'; %Posthoc correction after Kruskal Wallis non parametric anova

params.colors_audio_opto       = {[212 0 41] [255 122 94] [10 10 10]}; 
params.colors_audio_opto       = cellfun(@(x) x/256,params.colors_audio_opto,'UniformOutput',false);

params.colors_visual_opto       = {[40 0 150] [0 173 240] [10 10 10]}; 
params.colors_visual_opto       = cellfun(@(x) x/256,params.colors_visual_opto,'UniformOutput',false);

params.lines_audio_opto       = {'-o' '--o' ':o'}; 
params.lines_visual_opto       = {'-o' '--o' ':o'}; 

%% Filter sessions with optogenetic manipulation in V1 or PPC:
sesids            = unique(sessionData.session_ID(sessionData.UseOpto==1 & ...
    (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC')) & ...
    sessionData.Photostimpower >= params.minPhotoStimPower));
fprintf('Selected %d/%d sessions with optical manipulation in V1 or PPC\n',length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Trim post change delayed inhibition (for V1 - 2nd bump project):
trialfields     = fieldnames(trialData);
idx             = ~(trialData.PostChangeOptoStart > 0);
for iF = 1:length(trialfields)
    trialData.(trialfields{iF})    = trialData.(trialfields{iF})(idx);
end

%% Filter out sessions with only delayed or prechange photostimulation:
% trialData.PostChangeOptoStart(trialData.PostChangeOptoStart<0) = 0;

sesids              = unique(trialData.session_ID(trialData.PostChangeOptoStart==0));
fprintf('Removed %d/%d sessions with prechange or delayed opto inhibition\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Initialize structure for saving output fit parameters:
nSessions               = length(sessionData.session_ID);
dVis                    = NaN(nSessions,2,3); %Init matrix for storing all visual dprime data
dAud                    = NaN(nSessions,2,3); %Init matrix for storing all audio dprime data
cVis                    = NaN(nSessions,2,3); %Init matrix for storing all visual dprime data
cAud                    = NaN(nSessions,2,3); %Init matrix for storing all audio dprime data
TotalResp               = NaN(3,3,3,3,nSessions); %Init matrix for storing all hitrate data

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
    
    %Correction: If psychometric protocol, take intermediate values to compare with only 2 levels present
    if numel(visconditions) > 2
        FullRespMat = FullRespMat(:,[1 3 end],:,:);
        FullnTrialsMat = FullnTrialsMat(:,[1 3 end],:);
    end
    if numel(auconditions) > 2
        FullRespMat = FullRespMat([1 3 end],:,:,:);
        FullnTrialsMat = FullnTrialsMat([1 3 end],:,:);
    end
    
    FullRespMat(FullRespMat==0) = 1/50; %for dprime impossible to have zero or one, becomes infinite
    FullRespMat(FullRespMat==1) = 1 - 1/50;

    %Store values:
    TotalResp(:,:,:,:,iSes)     = FullRespMat;
    FullRespMat_trialn          = FullRespMat .* repmat(FullnTrialsMat,1,1,1,3);
    
    %Compute d-prime for each condition:
    for iTrial = 1:2
        for iOpto = 1:3
            outcome         = NaN(3,3);
            outcome(1,:)    = squeeze(FullRespMat_trialn(1,1+iTrial,iOpto,:)); %visual trials
            outcome(2,:)    = squeeze(FullRespMat_trialn(1+iTrial,1,iOpto,:)); %audio trials
            outcome(3,:)    = squeeze(FullRespMat_trialn(1,1,iOpto,:)); %probe trials
            outcome         = outcome(:,[2 1 3]); %Swap response coding, visual first.

            if nansum(outcome(:))>params.sumTrialsCondition
                [dVis(iSes,iTrial,iOpto),dAud(iSes,iTrial,iOpto),cVis(iSes,iTrial,iOpto),cAud(iSes,iTrial,iOpto)] = ...
                    MOL_Fit_2ADC_Full_Session(tempsessionData,trialData,0,outcome);
%                 dVis(iSes,iTrial,iOpto) = norminv(FullRespMat(1,1+iTrial,iOpto,2)) - norminv(FullRespMat(1,1,iOpto,2)); %quick approximation, but not accurate
%                 dAud(iSes,iTrial,iOpto) = norminv(FullRespMat(1+iTrial,1,iOpto,1)) - norminv(FullRespMat(1,1,iOpto,1)); %quick approximation, but not accurate
            end
        end
    end
end

fprintf('\n\n Finished fitting behavioral model.\n\n\n')

%% FIGURES:

%% Show average response rate figure for photoinhibition @V1 in MST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'V1') & ...
    ismember(sessionData.Experiment,{'BehaviorConflict' 'ChangeDetectionConflict'});

fprintf( '\nHit rates:\n')
MOL_plotOptoBehavior_Rates_PPC(params,TotalResp(:,:,:,:,idx_ses))
fprintf( '\nDprime:\n')
MOL_plotOptoBehavior_Dprime_PPC(params,dVis(idx_ses,:,:),dAud(idx_ses,:,:))

%% Show average response rate figure for photoinhibition @PPC in MST animals:
idx_ses = strcmp(sessionData.PhotostimArea,'PPC') & ...
    ismember(sessionData.Experiment,{'BehaviorConflict' 'ChangeDetectionConflict'});

fprintf( '\nHit rates:\n')
MOL_plotOptoBehavior_Rates_PPC(params,TotalResp(:,:,:,:,idx_ses))
fprintf( '\nDprime:\n')
MOL_plotOptoBehavior_Dprime_PPC(params,dVis(idx_ses,:,:),dAud(idx_ses,:,:))

%%



%% Additional analyses:  effect per mouse

uMice = unique(sessionData.mousename(strcmp(sessionData.PhotostimArea,'PPC'))); %get the 5 different mice with PPC inactivation
nMice = length(uMice);

dVis_mice = NaN(nMice,2,2); %dim2: thr/max: dim3: opto/no opto
dAud_mice = NaN(nMice,2,2); %dim2: thr/max: dim3: opto/no opto
for iM = 1:nMice %loop over mice:
    idx_ses = strcmp(sessionData.mousename,uMice{iM});
    dVis_mice(iM,:,:) = nanmean(dVis(idx_ses,:,1:2),1);
    dAud_mice(iM,:,:) = nanmean(dAud(idx_ses,:,1:2),1);
end

markers = {'h-' '^-' 'o-' 's-' 'd-'}; %different markers per mouse

%Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.5 .16 .3]); hold all;
handles = [];

for iM = 1:nMice
    temp = randn()*0.1; %Get random small horizontal offset
    %FOr visuaL:
    plot([1 2]+temp,[dVis_mice(iM,1,1) dVis_mice(iM,1,2)],'k-','LineWidth',1); %plot lines
    plot([1 2]+temp,[dVis_mice(iM,2,1) dVis_mice(iM,2,2)],'k-','LineWidth',1);
    
    plot(1+temp,dVis_mice(iM,1,1),markers{iM},'Color','k','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',1); %plot markers
    plot(2+temp,dVis_mice(iM,1,2),markers{iM},'Color','k','MarkerEdgeColor',params.colors_visual_opto{2},'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',1);
    
    handles(iM) = plot(1+temp,dVis_mice(iM,2,1),markers{iM},'Color','k','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',params.colors_visual_opto{1},'MarkerSize',10,'LineWidth',1);
    plot(2+temp,dVis_mice(iM,2,2),markers{iM},'Color','k','MarkerEdgeColor',params.colors_visual_opto{2},'MarkerFaceColor',params.colors_visual_opto{2},'MarkerSize',10,'LineWidth',1);
    %For auditory:
    plot([3 4]+temp,[dAud_mice(iM,1,1) dAud_mice(iM,1,2)],'k-','LineWidth',1);
    plot([3 4]+temp,[dAud_mice(iM,2,1) dAud_mice(iM,2,2)],'k-','LineWidth',1);
    
    plot(3+temp,dAud_mice(iM,1,1),markers{iM},'Color','k','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',1);
    plot(4+temp,dAud_mice(iM,1,2),markers{iM},'Color','k','MarkerEdgeColor',params.colors_audio_opto{2},'MarkerFaceColor',[1 1 1],'MarkerSize',10,'LineWidth',1);
    
    plot(3+temp,dAud_mice(iM,2,1),markers{iM},'Color','k','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',params.colors_audio_opto{1},'MarkerSize',10,'LineWidth',1);
    plot(4+temp,dAud_mice(iM,2,2),markers{iM},'Color','k','MarkerEdgeColor',params.colors_audio_opto{2},'MarkerFaceColor',params.colors_audio_opto{2},'MarkerSize',10,'LineWidth',1);
end
%Figure make up:
xlim([0.5 4.5])
ylim([-0.5 3])
ylabel('Dprime')
set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Ctrl' 'Opto' 'Ctrl' 'Opto'},'YTick',[0 1 2 3]);
legend(handles,uMice,'Location','NorthWest'); legend boxoff;

%% Additional analyses: bias effect per mouse
%Show relationship between bias/criterion and dprime effect:
figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.5 .13 .24]); hold all;

dAud_delta = dAud_mice(:,:,2) - dAud_mice(:,:,1);
dVis_delta = dVis_mice(:,:,2) - dVis_mice(:,:,1);

cVis_mice = NaN(nMice,2); 
cAud_mice = NaN(nMice,2); 
for iM = 1:nMice
    idx_ses = strcmp(sessionData.mousename,uMice{iM});
    cVis_mice(iM,:,:) = nanmean(cVis(idx_ses,:,1),1);
    cAud_mice(iM,:,:) = nanmean(cAud(idx_ses,:,1),1);
end

scatter(cVis_mice(:,1),dVis_delta(:,1),20,'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',[1 1 1]);
scatter(cVis_mice(:,2),dVis_delta(:,2),20,'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',params.colors_visual_opto{1});

scatter(cAud_mice(:,1),dAud_delta(:,1),20,'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',[1 1 1]);
scatter(cAud_mice(:,2),dAud_delta(:,2),20,'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',params.colors_audio_opto{1});
%Figure makeuP:
xlabel('Criterion')
ylabel('Delta d-prime')
xlim([-1 2])
ylim([-2.2 2.2])

%% Additional analyses: bias effect per session
idx_ses = strcmp(sessionData.PhotostimArea,'PPC') & ...
    ismember(sessionData.Experiment,{'BehaviorConflict' 'ChangeDetectionConflict'});

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.5 .14 .21]); hold all;
dAud_delta = dAud(:,:,2) - dAud(:,:,1);
dVis_delta = dVis(:,:,2) - dVis(:,:,1);

markers = {'h' '^' 'o' 's' 'd'}; %different markers per mouse

for iM = 1:nMice
    idx_all = idx_ses & strcmp(sessionData.mousename,uMice{iM});
    scatter(cVis(idx_all,1,1),dVis_delta(idx_all,1),20,markers{iM},'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',[1 1 1]);
    scatter(cVis(idx_all,2,1),dVis_delta(idx_all,2),20,markers{iM},'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',params.colors_visual_opto{1});
end

% scatter(cVis(idx_ses,1,1),dVis_delta(idx_ses,1),20,'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',[1 1 1]);
% scatter(cVis(idx_ses,2,1),dVis_delta(idx_ses,2),20,'filled','MarkerEdgeColor',params.colors_visual_opto{1},'MarkerFaceColor',params.colors_visual_opto{1});
xdata = squeeze(cVis(idx_ses,1,1:2)); 
ydata = dVis_delta(idx_ses,1:2);
idx = ~isnan(xdata) & ~isnan(ydata);
xdata = xdata(idx);
ydata = ydata(idx);

[bf10,r,~] = bf.corr(xdata(:),ydata(:));
fprintf('\nWe found no significant correlation between criterion\n and delta dprime (Bayesian correlation; visual: r=%1.2f, BF=%1.2f; ',r,bf10)
h = bayesstar([0.5 1],bf10);
text(0.75,1,sprintf('r=%1.2f',r),'Color',[0 0 0.8])
set(h(2),'Color',[0 0 0.8])

for iM = 1:nMice
    idx_all = idx_ses & strcmp(sessionData.mousename,uMice{iM});
    scatter(cAud(idx_all,1,1),dAud_delta(idx_all,1),20,markers{iM},'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',[1 1 1]);
    scatter(cAud(idx_all,2,1),dAud_delta(idx_all,2),20,markers{iM},'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',params.colors_audio_opto{1});
end
% 
% scatter(cAud(idx_ses,1,1),dAud_delta(idx_ses,1),20,'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',[1 1 1]);
% scatter(cAud(idx_ses,2,1),dAud_delta(idx_ses,2),20,'filled','MarkerEdgeColor',params.colors_audio_opto{1},'MarkerFaceColor',params.colors_audio_opto{1});
xdata = squeeze(cAud(idx_ses,1,1:2)); 
ydata = dAud_delta(idx_ses,1:2);
idx = ~isnan(xdata) & ~isnan(ydata);
xdata = xdata(idx);
ydata = ydata(idx);

[bf10,r,~] = bf.corr(xdata(:),ydata(:));
fprintf('auditory: r=%1.2f, BF=%1.2f)\n',r,bf10)
text(0.95,0.8,sprintf('r=%1.2f',r),'Color',[0.8 0 0])
h = bayesstar([1 1.5],bf10);
set(h(2),'Color',[0.8 0 0])

%Figure makeup:
xlabel('Criterion')
ylabel('Delta d-prime')
xlim([-1 2])
ylim([-2.2 2.2])
ylim([-1.5 1.5])

%% Additional analyses: reaction time effects

trial_sesidx            =  ismember(trialData.session_ID,sessionData.session_ID(sessionData.UseOpto==1 & strcmp(sessionData.PhotostimArea,'PPC') & sessionData.Photostimpower >= params.minPhotoStimPower));

datatoplot = NaN(10000,4);
idx = trial_sesidx & trialData.vecResponse==1 & ~(trialData.hasphotostim==1);
datatoplot(1:sum(idx),1)       = trialData.responseLatency(idx);

idx = trial_sesidx & trialData.vecResponse==1 & trialData.hasphotostim==1;
datatoplot(1:sum(idx),2)       = trialData.responseLatency(idx);

idx = trial_sesidx & trialData.vecResponse==2 & ~(trialData.hasphotostim==1);
datatoplot(1:sum(idx),3)       = trialData.responseLatency(idx);

idx = trial_sesidx & trialData.vecResponse==2 & trialData.hasphotostim==1;
datatoplot(1:sum(idx),4)       = trialData.responseLatency(idx);

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .17 .47]); hold all;
boxplot(datatoplot, 'plotstyle','compact','Colors',[0.7 0 0; 1 0.3 0.3; 0 0 0.7; 0.3 0.3 1])

%Classical statistics (NHST):
% p = ranksum(datatoplot(:,1),datatoplot(:,2));
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions
% text(1.2,200e3,sprintf('p=%1.3f',p),'FontSize',15)
% 
% p = ranksum(datatoplot(:,3),datatoplot(:,4));
% sigstar([3 4],p) %use sigstar function to identify signficantly different conditions
% text(3.2,200e3,sprintf('p=%1.3f',p),'FontSize',15)

%Bayesian statitics:
tempbf      = bf.ttest2(datatoplot(:,1),datatoplot(:,2));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([1 2]),nanmean(datatoplot(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(datatoplot(:,1),datatoplot(:,2),'paired');
fprintf('Visual hits: %3.0f ms (control), %3.0f ms (opto), %d vs %d trials, Cohen''s d = %1.3f, BF=%3.3f\n',nanmedian(datatoplot(:,1))*1e-3,nanmedian(datatoplot(:,2))*1e-3,sum(~isnan(datatoplot(:,1))),sum(~isnan(datatoplot(:,2))),tempd,tempbf)

tempbf      = bf.ttest2(datatoplot(:,3),datatoplot(:,4));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([3 4]),nanmean(datatoplot(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(datatoplot(:,3),datatoplot(:,4),'paired');
fprintf('Audio hits: %3.0f ms (control), %3.0f ms (opto), %d vs %d trials, Cohen''s d = %1.3f, BF=%3.3f\n',nanmedian(datatoplot(:,3))*1e-3,nanmedian(datatoplot(:,4))*1e-3,sum(~isnan(datatoplot(:,3))),sum(~isnan(datatoplot(:,4))),tempd,tempbf)

%Make up:
set(gca, 'XTick', 1:4,'XTickLabels', {'Control' 'Photostim' 'Control' 'Photostim'},'XTickLabelRotation',45)
set(gca, 'YTick', (0:250:1200)*1e3,'YTickLabels',0:250:1200)
ylim([0e3 1200e3])
xlim([0.7 4.3])
ylabel('Reaction Time')

%% Additional analysis: effect on conflict trials:

trial_sesidx            =  ismember(trialData.session_ID,sessionData.session_ID(sessionData.UseOpto==1 & strcmp(sessionData.PhotostimArea,'PPC') & sessionData.Photostimpower >= params.minPhotoStimPower));

datatoplot = NaN(10000,4);
idx = trial_sesidx & strcmp(trialData.trialType,'C') & trialData.vecResponse==1 & ~(trialData.hasphotostim==1);
datatoplot(1:sum(idx),1)       = trialData.responseLatency(idx);

idx = trial_sesidx & strcmp(trialData.trialType,'C') & trialData.vecResponse==1 & trialData.hasphotostim==1;
datatoplot(1:sum(idx),2)       = trialData.responseLatency(idx);

idx = trial_sesidx & strcmp(trialData.trialType,'C') & trialData.vecResponse==2 & ~(trialData.hasphotostim==1);
datatoplot(1:sum(idx),3)       = trialData.responseLatency(idx);

idx = trial_sesidx & strcmp(trialData.trialType,'C') & trialData.vecResponse==2 & trialData.hasphotostim==1;
datatoplot(1:sum(idx),4)       = trialData.responseLatency(idx);

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .17 .47]); hold all;
title('Conflict trials')
boxplot(datatoplot, 'plotstyle','compact','Colors',[0.7 0 0; 1 0.3 0.3; 0 0 0.7; 0.3 0.3 1])

%Classical statistics (NHST):
% p = ranksum(datatoplot(:,1),datatoplot(:,2));
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions
% text(1.2,200e3,sprintf('p=%1.3f',p),'FontSize',15)
% 
% p = ranksum(datatoplot(:,3),datatoplot(:,4));
% sigstar([3 4],p) %use sigstar function to identify signficantly different conditions
% text(3.2,200e3,sprintf('p=%1.3f',p),'FontSize',15)

%Bayesian statitics:
tempbf      = bf.ttest2(datatoplot(:,1),datatoplot(:,2));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([1 2]),nanmean(datatoplot(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(datatoplot(:,1),datatoplot(:,2),'paired');
% fprintf('Bayesian ttest (RT conflict audio choice): %d ctrl trials, %d opto trials, Cohen''s d = %1.3f, BF10=%3.3f\n',sum(~isnan(datatoplot(:,1))),sum(~isnan(datatoplot(:,2))),tempd,tempbf)
fprintf('Auditory choice trials: %3.0f ms (control), %3.0f ms (opto), %d vs %d trials, Cohen''s d = %1.3f, BF=%3.3f\n',nanmedian(datatoplot(:,1))*1e-3,nanmedian(datatoplot(:,2))*1e-3,sum(~isnan(datatoplot(:,1))),sum(~isnan(datatoplot(:,2))),tempd,tempbf)

tempbf      = bf.ttest2(datatoplot(:,3),datatoplot(:,4));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([3 4]),nanmean(datatoplot(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(datatoplot(:,3),datatoplot(:,4),'paired');
% fprintf('Bayesian ttest (RT conflict visual choice): %d ctrl trials, %d opto trials, Cohen''s d = %1.3f, BF10=%3.3f\n',sum(~isnan(datatoplot(:,3))),sum(~isnan(datatoplot(:,4))),tempd,tempbf)
fprintf('Visual choice trials: %3.0f ms (control), %3.0f ms (opto), %d vs %d trials, Cohen''s d = %1.3f, BF=%3.3f\n',nanmedian(datatoplot(:,3))*1e-3,nanmedian(datatoplot(:,4))*1e-3,sum(~isnan(datatoplot(:,3))),sum(~isnan(datatoplot(:,4))),tempd,tempbf)

%Make up:
set(gca, 'XTick', 1:4,'XTickLabels', {'Control' 'Photostim' 'Control' 'Photostim'},'XTickLabelRotation',45)
set(gca, 'YTick', (0:250:1200)*1e3,'YTickLabels',0:250:1200)
ylim([0e3 1200e3])
xlim([0.7 4.3])
ylabel('Reaction Time')

%% Additional analysis: effect on conflict trials:

sesselec            =  sessionData.session_ID(sessionData.UseOpto==1 & strcmp(sessionData.PhotostimArea,'PPC') & sessionData.Photostimpower >= params.minPhotoStimPower);

datatoplot = NaN(3,3,length(sesselec));

for ses = 1:length(sesselec)
    sesid                                                       = sesselec(ses);
    [tempsessionData,temptrialData]                             = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually:
    [visconditions,auconditions,FullRespMat,FullnTrialsMat]     = MOL_Psy_GetTrialMatOpto(tempsessionData,temptrialData);
    datatoplot(:,:,ses)     = squeeze(nansum(nansum(FullRespMat(2:end,2:end,:,:),2),1));
end

ratioA = squeeze(datatoplot(1:2,1,:)) ./ squeeze(sum(datatoplot(1:2,:,:),2));

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.1 .17 .47]); hold all;
for ses = 1:length(sesselec)
    plot([1 2],[ratioA(1,ses) ratioA(2,ses)],'.-','Color','k','MarkerSize',30,'LineWidth',1); hold all;
end

ratioA = ratioA';

%Classical statitics:
% p = ranksum(ratio(:,1),ratio(:,2));
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions
% text(1.2,0.1,sprintf('p=%1.3f',p),'FontSize',15)

%Bayesian statitics:
tempbf      = bf.ttest(ratioA(:,1),ratioA(:,2));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([1 2]),nanmean(ratioA(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(ratioA(:,1),ratioA(:,2),'paired');
fprintf('Bayesian ttest (Auditory choice conflict trials): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',size(ratioA,1),tempd,tempbf)

%Make up:
plot([-1 5],[0.5 0.5],'k:','LineWidth',0.5)
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0 1])
xlim([0.7 2.3])
ylabel('%Choice auditory')

ratioV = squeeze(datatoplot(1:2,2,:)) ./ squeeze(sum(datatoplot(1:2,:,:),2));

figure; set(gcf,'color','w','units','normalized','Position', [0.7 0.1 .17 .47]); hold all;
for ses = 1:length(sesselec)
    plot([1 2],[ratioV(1,ses) ratioV(2,ses)],'.-','Color','k','MarkerSize',30,'LineWidth',1); hold all;
end
ratioV = ratioV';

%Classical statitics:
% p = ranksum(ratio(:,1),ratio(:,2));
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions
% text(1.2,0.1,sprintf('p=%1.3f',p),'FontSize',15)

%Bayesian statitics:
tempbf      = bf.ttest(ratioV(:,1),ratioV(:,2));
bfsymb      = MOL_BFtoSymbol(tempbf);
text(mean([1 2]),nanmean(ratioV(:)),bfsymb,'FontSize',15)
tempd       = computeCohen_d(ratioV(:,1),ratioV(:,2),'paired');
fprintf('Bayesian ttest (Visual choice conflict trials): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',size(ratioV,1),tempd,tempbf)

%Make up:
plot([-1 5],[0.5 0.5],'k:','LineWidth',0.5)
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0 1])
xlim([0.7 2.3])
ylabel('%Choice Visual')


