%% MOL_GLM_History

%% Get Data: 
animals = {'2003'    '2004'    '2007'    '2009'    '2010'    '2011'   '2012'    '2013'    '2019'    '2020'    '2021'    '2022'    '2023'    '2026'    '2027'    '2030'    '2031'};
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict' 'ChangeDetectionConflict'},animals,[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Filter out sessions with muscimol:
sesids              = sessionData.session_ID(cellfun(@isempty,sessionData.MuscimolArea));
fprintf('Removed %d/%d sessions with muscimol manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% Filter out sessions with photostimulation in V1 or PPC:
sesids              = sessionData.session_ID(~(sessionData.UseOpto & (strcmp(sessionData.PhotostimArea,'V1') | strcmp(sessionData.PhotostimArea,'PPC'))));
fprintf('Removed %d/%d sessions with optogenetic manipulations\n',length(sessionData.session_ID)-length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%% 

params.kfold            = 3;

% params.glmmethod        = 'mnrfit';
params.glmmethod        = 'glmnet';

params.lambdastring     = 'lambda_min'; %'lambda_1se'
params.elastic_alpha    = 0.95;

% params.posthoctest      = 'tukey-kramer';

%% Loop over mice and then sessions:

mouseids                = unique(sessionData.mousename)';
nMice                   = length(mouseids);
B                       = {};
CV_Perf                 = NaN(10,length(mouseids));
% Y_Hat                   = NaN(10,length(mouseids),5000);

fprintf('Fitting data on ~%2.2f concatenated sessions per mouse\n',length(sessionData.session_ID)/nMice)

for iMou = 1:nMice
    fprintf('Fitting GLM behavioral models for mouse %d/%d\n',iMou,length(mouseids))
    
    %Select all the data from this mouse that fits the criteria:
    sesselec                = sessionData.session_ID(strcmp(sessionData.mousename,mouseids(iMou)));
%     [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    
    X_rand              = rand(size(temptrialData.stimChange));
    
    % Session duration = trial number:
    X_trialNumber       = temptrialData.trialNum;
    
    %Modality (trial type):
    X_trialType  = NaN(size(temptrialData.session_ID));
    X_trialType(strcmp(temptrialData.trialType,'Y')) = 1;
    X_trialType(strcmp(temptrialData.trialType,'X')) = 2;
    X_trialType(strcmp(temptrialData.trialType,'P')) = 3;
    X_trialType(strcmp(temptrialData.trialType,'C')) = 4;
    if any(isnan(X_trialType))
        error('unknown trial types')
    end
    
    %Visual Intensity:
    X_visualOriChange       = zeros(size(temptrialData.session_ID));
    X_visualOriChange(strcmp(temptrialData.trialType,'X')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'X')));
    X_visualOriChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'C')));
    %Log Visual Intensity:
    X_logvisualOriChange    = X_visualOriChange;
    X_logvisualOriChange(X_logvisualOriChange==0) = 0.1;
    X_logvisualOriChange    = log(X_logvisualOriChange);
    
    %Auditory Intensity:
    X_audioFreqChange       = zeros(size(temptrialData.session_ID));
    X_audioFreqChange(strcmp(temptrialData.trialType,'Y')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'Y')));
    X_audioFreqChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'C')));
    %Log Auditory Intensity:
    X_logaudioFreqChange    = X_audioFreqChange;
    X_logaudioFreqChange(X_logaudioFreqChange==0) = 0.1;
    X_logaudioFreqChange    = log(X_logaudioFreqChange);
    
    %Rewards: (binary)
    X_Correct       = temptrialData.correctResponse;
    if tempsessionData.VisualLeftCorrectSide
        X_visualCorrect     = double(temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L'));
        X_audioCorrect  	= double(temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R'));
    else
        X_visualCorrect     = double(temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R'));
        X_audioCorrect      = double(temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L'));
    end
    
    X_LastLick      = NaN(size(X_Correct));
    for iTrial = 1:length(temptrialData.lickSide)
        templick                    = [temptrialData.lickSide{max([1 iTrial-3]):iTrial}];
        templick                    = templick([temptrialData.lickTime{max([1 iTrial-3]):iTrial}] < temptrialData.stimChange(iTrial));
        if ~isempty(templick)
            X_LastLick(iTrial)          = templick(end);
        end
    end
    if tempsessionData.VisualLeftCorrectSide
        X_LastLick(X_LastLick==76) = 0;     X_LastLick(X_LastLick==82) = 1; %assign left and right licks different values
    else
        X_LastLick(X_LastLick==76) = 1;     X_LastLick(X_LastLick==82) = 0; %assign left and right licks different values
    end
    X_LastLick(isnan(X_LastLick))       = randi(2);
    
    %Photostim:
    if isfield(temptrialData,'hasphotostim')
        X_photoStim         = temptrialData.hasphotostim;
    else X_photoStim = rand(size(X_Correct));
    end
    
    %Response vector:
    %Auditory response  = 1
    %Visual response    = 2
    %No response        = 3
    Y_Response = temptrialData.vecResponse;

    iModel = 0;
    
    %Null model: random variable
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Null'; %#ok<*SAGROW>
    X                           = [X_rand]; %#ok<*NBRAK>
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    %Session duration model: trial number
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Trial';
    X                           = [X_trialNumber];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    %Sensory history model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory history';
    X                           = [[0; X_logvisualOriChange(1:end-1)] [0; X_logaudioFreqChange(1:end-1)]];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);

    % Reward History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Reward history';
    X                           = [[0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)]];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % Motor history model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Motor history';
    X                           = [X_LastLick];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    %Log Sensory model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory';
    X                           = [X_logvisualOriChange X_logaudioFreqChange];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % Sensory + History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory + Motor History';
    X                           = [X_trialNumber X_logvisualOriChange X_logaudioFreqChange X_LastLick];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % Sensory Motor History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Full';
    X                           = [X_trialNumber X_LastLick [0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)] [0; X_logvisualOriChange(1:end-1)] [0; X_logaudioFreqChange(1:end-1)] X_logvisualOriChange X_logaudioFreqChange];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);

    %     Export data to .csv file
    %     X                   = [X_visualOriChange X_audioFreqChange X_visualCorrect X_audioCorrect];
    %     rootdir             = 'E:\Documents\PhD\TempExportPietro\';
    %     csvwrite(sprintf('%s%s.csv',rootdir,mouseids{iMou}),[Y_Response X],0) %writes matrix M to file filename as comma-separated values.
    
end

nModels = size(B,1); %count number of different models that were fit
CV_Perf = CV_Perf(1:nModels,:); %trim the cross validated performance

%% 

legendhandles = [];
% Make figure of cross-validated performance on test data:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.35 .27 .55]); hold all;
xlim([0.5 nModels+0.5])

% Colors  = parula(length(mouseids));
% Colors  = lines(length(mouseids));

% cmp = getPyPlot_cMap('GnBu_r',128,keepAlpha,pyCmd)
Colors = getPyPlot_cMap('brg',length(mouseids));
Colors = getPyPlot_cMap('brg',5);

xloc    = 1; %position of performance per mouse/mean

for iMod=1:nModels
    for iMou = 1:nMice
        h = plot(xloc-0.4+rand()*0.2,CV_Perf(iMod,iMou),'.','Color',Colors(iMou,:),'MarkerSize',40); hold all;
        if iMod==1
            legendhandles(iMou) = h;
        end
    end
    mPerf        = nanmean(CV_Perf(iMod,:));
    semPerf      = nanstd(CV_Perf(iMod,:)) / sqrt(size(CV_Perf,2));
    h = errorbar(xloc,mPerf,semPerf,semPerf,'k','LineWidth',2);
%     h.CapSize = 12;
    plot(xloc,mPerf,'k.','MarkerSize',30)
    xloc = xloc+1;
end

%Statistical testing:
%perform Friedman nonparametric repeated measures anova:
% [~,~,stats]         = friedman(CV_Perf',1,'off');
% comptable           = multcompare(stats,'display','off','alpha',0.05,'ctype',params.posthoctest); %do posthoc
% comptable           = comptable(comptable(:,end)<0.05,:); %Filter only significant
% sigstar(mat2cell(comptable(:,1:2),ones(size(comptable,1),1)),comptable(:,end)) %use sigstar function to identify

%Statistical testing paired testing, pseudobonferroni corrected:
for iM = 1:nModels-1
    p = signrank(CV_Perf(iM,:),CV_Perf(iM+1,:));
    if p*nModels < 0.05
        sigstar([iM iM+1],p*nModels) %use sigstar function to identify
    end
end

%Make up:
set(gca, 'XTick', 1:nModels,'XTickLabels', Model_string,'XTickLabelRotation',45)
ylim([0.33 0.75])
xlim([0.2 nModels+0.2])
ylabel('Cross-validated % Correct')
legend(legendhandles,mouseids,'Location','NorthWest','FontSize',8); legend boxoff;

variance = reshape(CV_Perf',nMice*nModels,1);
groups = reshape(repmat(1:nModels,nMice,1),nMice*nModels,1);
tbl = table(variance,groups);

tempBF = bf.anova(tbl,'variance~groups');

%%
for iMod = 2:6
    tempbf = bf.ttest(CV_Perf(1,:),CV_Perf(iMod,:));    
    tempd = computeCohen_d(CV_Perf(iMod+1,:),CV_Perf(1,:),'paired');
    fprintf('%s: Cohens'' d: %2.2f,BF=%3.2f; ',Model_string{iMod},tempd,tempbf);
end
fprintf('\n\n' )

%% 
nModels = size(B,1);
legendhandles = [];
% Make figure of cross-validated performance on test data:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .8 .8]); hold all;

Colors  = parula(length(mouseids));
xloc    = 1; %position of performance per mouse/mean

for iMod=1:nModels
    for iMou = 1:length(mouseids)
        h = plot(xloc,CV_Perf(iMod,iMou),'.','Color',Colors(iMou,:),'MarkerSize',40); hold all;
        xloc = xloc+1;
    end
    
    mPerf        = nanmean(CV_Perf(iMod,:));
    semPerf      = nanstd(CV_Perf(iMod,:)) / sqrt(size(CV_Perf,2));
    plot(xloc,mPerf,'k.','MarkerSize',30)
    legendhandles(iMod) = errorbar(xloc,mPerf,semPerf,semPerf,'k','LineWidth',2);
    xloc = xloc+2;
end

%Make up:
set(gca, 'XTick', 1:(length(mouseids)+2)*nModels,'XTickLabels', repmat([mouseids 'Mean' ' '],1,nModels),'XTickLabelRotation',45)
ylim([0.33 0.7])
xlim([0.5 (length(mouseids)+2)*nModels+0.5])
ylabel('Cross-validated % Correct')
legend(legendhandles,Model_string,'Location','NorthWest');


%% 
% Make example of the fitting: 


params.examplemouse    = '2003';

%Select all the data from this mouse that fits the criteria:
sesselec                = sessionData.session_ID(strcmp(sessionData.mousename,params.examplemouse));
[tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:

% Session duration = trial number:
X_trialNumber       = temptrialData.trialNum;

%Visual Intensity:
X_visualOriChange       = zeros(size(temptrialData.session_ID));
X_visualOriChange(strcmp(temptrialData.trialType,'X')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'X')));
X_visualOriChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.visualOriChange(strcmp(temptrialData.trialType,'C')));
%Log Visual Intensity:
X_logvisualOriChange    = X_visualOriChange;
X_logvisualOriChange(X_logvisualOriChange==0) = 0.1;
X_logvisualOriChange    = log(X_logvisualOriChange);

%Auditory Intensity:
X_audioFreqChange       = zeros(size(temptrialData.session_ID));
X_audioFreqChange(strcmp(temptrialData.trialType,'Y')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'Y')));
X_audioFreqChange(strcmp(temptrialData.trialType,'C')) = abs(temptrialData.audioFreqChange(strcmp(temptrialData.trialType,'C')));
%Log Auditory Intensity:
X_logaudioFreqChange    = X_audioFreqChange;
X_logaudioFreqChange(X_logaudioFreqChange==0) = 0.1;
X_logaudioFreqChange    = log(X_logaudioFreqChange);

%Rewards: (binary)
X_Correct       = temptrialData.correctResponse;
if tempsessionData.VisualLeftCorrectSide
    X_visualCorrect     = double(temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L'));
    X_audioCorrect  	= double(temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R'));
else
    X_visualCorrect     = double(temptrialData.rightCorrect & strcmp(temptrialData.responseSide,'R'));
    X_audioCorrect      = double(temptrialData.leftCorrect & strcmp(temptrialData.responseSide,'L'));
end

X_LastLick      = NaN(size(X_Correct));
for iTrial = 1:length(temptrialData.lickSide)
    templick                    = [temptrialData.lickSide{1:iTrial}];
    templick                    = templick([temptrialData.lickTime{1:iTrial}] < temptrialData.stimChange(iTrial));
    if ~isempty(templick)
        X_LastLick(iTrial)          = templick(end);
    end
end
X_LastLick(isnan(X_LastLick))       = randi(2);
X_LastLick(X_LastLick==76) = 1;     X_LastLick(X_LastLick==82) = 2; %assign left and right licks different values

%Response vector:
%Auditory response  = 1
%Visual response    = 2
%No response        = 3
Y_Response = temptrialData.vecResponse;

% Sensory Motor History model:
X                           = [X_trialNumber X_LastLick [0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)] X_logvisualOriChange X_logaudioFreqChange];
[~,~,Y_hat_choice]     = MOL_fitGLM_session(X,Y_Response,params);

% X                           = [rand(size(X_visualOriChange))]; %#ok<*NBRAK>
% [~,~,Y_hat_choice_rand]     = MOL_fitGLM_session(X,Y_Response,params);

%% 

%Swap values: 
Y_Response(Y_Response==3) = 4; Y_Response(Y_Response==2) = 3; 
Y_Response(Y_Response==4) = 2; 

Y_hat_choice(Y_hat_choice==3) = 4; Y_hat_choice(Y_hat_choice==2) = 3; 
Y_hat_choice(Y_hat_choice==4) = 2; 
% 
% Y_hat_choice_rand(Y_hat_choice_rand==3) = 4; Y_hat_choice_rand(Y_hat_choice_rand==2) = 3; 
% Y_hat_choice_rand(Y_hat_choice_rand==4) = 2; 

%% 
chunksize = 30;

fitqual     = Y_Response==Y_hat_choice; 
nParts      = length(fitqual)-mod(length(fitqual),chunksize);

fitqual     = reshape(fitqual(1:nParts),chunksize,nParts/chunksize);
fitqual     = sum(fitqual,1);
[~,idx]     = max(fitqual);
params.exampletrials = idx*chunksize:idx*chunksize+chunksize;

params.exampletrials  = params.exampletrials  - chunksize;

figure; set(gcf,'color','w','units','normalized','Position', [0.5 0.45 0.4 0.35]); hold all;
plot(params.exampletrials+0.2,Y_hat_choice(params.exampletrials),'Color',[0.1 0.6 0.1],'LineWidth',3);
plot(params.exampletrials,Y_Response(params.exampletrials),'Color',[0 0 0],'LineWidth',3);

plot(params.exampletrials+0.2,Y_hat_choice(params.exampletrials),'.','Color',[0.1 0.8 0.1],'MarkerSize',30);
plot(params.exampletrials,Y_Response(params.exampletrials),'k.','MarkerSize',30);

set(gca,'YTick',[1 2 3],'YTickLabels',{'Audio spout' 'No lick' 'Visual spout'})
xlim([params.exampletrials(1) params.exampletrials(end)])
ylim([0.8 3.2])

legend('Model','Data')

%%








%% Make figure of coefficients:
[~,nMice,~]  = size(B);
Modality_string = {'Auditory' 'Visual'};

iMod = 8;
nBeta       = size(B{iMod,1,:},1)-1;

% Predictors_string = {'TrialNum' 'MotorHist' 'VisHitHist' 'AudioHitHist' 'Visual Stim' 'Audio Stim'};
Predictors_string = {'TrialNum' 'MotorHist' 'VisHitHist' 'AudioHitHist' 'Vis Hist' 'Audio Hist' 'Visual Stim' 'Audio Stim'};

figure; set(gcf,'color','w','units','normalized','Position', [0.025 0.2 0.2*nBeta 0.65]); hold all;
for iModality = 1:2
    for iB = 1:nBeta
%         subplot(2,nBeta,iModality+(iB-1)*2)
        subplot(2,nBeta,iB+(iModality-1)*nBeta)
        for iMou = 1:nMice
            allmodelcoef            = [B{iMod,iMou}];
            nSes = size(allmodelcoef,2)/2;
            plot(allmodelcoef(iB+1,iModality:2:nSes*2),repmat(iMou,1,nSes)/10,'.','Color',Colors(iMou,:),'MarkerSize',45); hold all;
        end
        set(gca, 'YTick', [1:nMice]/10,'YTickLabels', mouseids,'YTickLabelRotation',15,'Fontsize',10)
        ylim([0 (nMice+1)/10])
        title(sprintf('%s on %s',Predictors_string{iB},Modality_string{iModality}),'Fontsize',10)
        plot([0,0],[0 (nMice+1)/10],'k:')
        temp = get(gca,'XLim');
        xlim([min(temp(1),-1) max(temp(2),1)])
    end
end
tightfig();


%% Make figure of coefficients:
[~,nMice,~]  = size(B);
Modality_string = {'Auditory' 'Visual'};

iModel = 8;
nBeta       = size(B{iModel,1,:},1)-1;

Predictors_string = {'TrialNum' 'Choice-1' 'VisHit-1' 'AudioHit-1' 'Visual Stim -1' 'Audio Stim -1' 'Visual Stim' 'Audio Stim'};

figure; set(gcf,'color','w','units','normalized','Position', [0.25 0.2 0.45 0.2]); hold all;
for iModality = 1:2
    subplot(1,2,iModality); hold all;
    title(sprintf('Weight on %s lick',Modality_string{iModality}))
    datatoplot            = cat(3,B{iModel,:});
    meantoplot            = nanmean(datatoplot(2:end,iModality,:),3);
%     errortoplot           = nanstd(datatoplot(2:end,iModality,:),[],3) / sqrt(nMice);
    errortoplot           = nanstd(datatoplot(2:end,iModality,:),[],3);
    
    errorbarxy(meantoplot,transpose(1:nBeta),errortoplot,errortoplot,zeros(nBeta,1),zeros(nBeta,1),{'k.', 'k', 'k'})
    plot(meantoplot,1:length(meantoplot),'.','Color',[0 0 0],'MarkerSize',25); hold all;
    
     for iVar = 1:nBeta
        tempBF = bf.ttest(datatoplot(iVar+1,iModality,:),0);
        text(-2,iVar,MOL_BFtoSymbol(tempBF))
    end
    plot([0,0],[0 10],'k:')

    ylim([0.5 8.5])
    xlim([-4 4])
    if iModality == 1
        set(gca,'YTick',1:8,'YTickLabels',Predictors_string)
    else set(gca,'YTick',[])
    end
    
end
tightfig();

