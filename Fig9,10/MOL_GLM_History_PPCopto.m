%% MOL_GLM_History
% PPC inactivation

%% Run GLM for PPC inhibition:
startover

%% Get Data: 
[Data] = MOL_GetData('E:','CHDET',{'BehaviorConflict' 'ChangeDetectionConflict'},{'2009' '2010' '2011' '2012' '2013'},[],{'sessionData' 'trialData'});
sessionData     = Data.sessionData;
trialData       = Data.trialData;

%% Remove last 20 trials:
trialData = MOL_RemoveLastnTrials(trialData,20);

%% Select sessions with photostimulation in PPC: 
sesids              = sessionData.session_ID(sessionData.UseOpto & strcmp(sessionData.PhotostimArea,'PPC') & sessionData.Photostimpower >= 2);
fprintf('Included %d/%d sessions with Opto in PPC\n',length(sesids),length(sessionData.session_ID));
[sessionData, trialData]        = MOL_getTempPerSes(sesids,sessionData,trialData);

%%
% params.removeConflict   = 1;

params.kfold            = 3;

params.glmmethod        = 'glmnet';

params.lambdastring     = 'lambda_min'; %'lambda_1se'
params.elastic_alpha    = 0.95;

%% Loop over mice:

uMice                   = unique(sessionData.mousename)';
nMice                   = length(uMice);
B                       = {};
CV_Perf                 = NaN(10,nMice);
% Y_Hat                   = NaN(10,nMice,5000);

for iMou = 1:nMice
    fprintf('Fitting GLM behavioral models for mouse %d/%d\n',iMou,nMice)
    
    %Select all the data from this mouse that fits the criteria:
    sesselec                = sessionData.session_ID(strcmp(sessionData.mousename,uMice(iMou)));
%     [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    [tempsessionData,temptrialData] = MOL_getTempPerSes(sesselec,sessionData,trialData);%Get the sessionData for each mouse individually:
    
    %Random predictor:
    X_rand              = rand(size(temptrialData.session_ID));

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
        templick                    = [temptrialData.lickSide{1:iTrial}];
        templick                    = templick([temptrialData.lickTime{1:iTrial}] < temptrialData.stimChange(iTrial));
        if ~isempty(templick)
            X_LastLick(iTrial)          = templick(end);
        end
    end
    X_LastLick(isnan(X_LastLick))       = randi(2);
    X_LastLick(X_LastLick==76) = 1;     X_LastLick(X_LastLick==82) = 2; %assign left and right licks different values
    
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
    
    % Full model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Full';
    X                           = [X_rand X_trialNumber X_LastLick [0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)] X_logvisualOriChange X_logaudioFreqChange];
    
    allX = X;
    Y = Y_Response;
    
    K                   = numel(unique(Y));             %Number of outcome levels (categorical response options)
    [N,P]               = size(allX);                   %N = num observations, P = num predictors
    
    %Normalize all predictors:
    allX = allX - repmat(min(allX,[],1),N,1);
    allX = allX ./ repmat(max(allX,[],1),N,1);
    
    allX(isnan(allX)) = 0;
    
    switch params.glmmethod
        case 'mnrfit'
            %not implemented
            
        case 'glmnet'
            %Get correct regularization parameter with crossval:
            glmoptions              = struct();
            glmoptions.alpha        = params.elastic_alpha;                  %(elastic-net mixing parameter)
            cvfit                   = cvglmnet(allX,Y,'multinomial',glmoptions,'deviance',params.kfold,[],false);
            glmoptions.lambda       = cvfit.lambda_min;
            
            %Analysis 1: fit on opto trials and fit on control trials and compare weights to see differences: (Akrami et al. 2018)
            trialidx_opto           = find(X_photoStim);
            temp                    = find(~X_photoStim);
            trialidx_ctrl           = temp(randperm(length(temp),length(trialidx_opto))); %take equal amount of trials in control condition
            %Actual fit without crossvalidation for control trials:
            glmfit                  = glmnet(allX(trialidx_ctrl,:),Y(trialidx_ctrl),'multinomial',glmoptions);
            Y_hat_choice            = glmnetPredict(glmfit,allX(trialidx_ctrl,:),glmoptions.lambda,'class');
            CV_ctrl_1(iMou)         = sum(Y_hat_choice == Y(trialidx_ctrl)) / numel(Y(trialidx_ctrl)); %#ok<*SAGROW>
            B_ctrl_1{iMou}          = cell2mat(glmfit.beta);
            %Actual fit without crossvalidation for opto trials:
            glmfit                  = glmnet(allX(trialidx_opto,:),Y(trialidx_opto),'multinomial',glmoptions);
            Y_hat_choice            = glmnetPredict(glmfit,allX(trialidx_opto,:),glmoptions.lambda,'class');
            CV_opto_1(iMou)         = sum(Y_hat_choice == Y(trialidx_opto)) / numel(Y(trialidx_opto));
            B_opto_1{iMou}          = cell2mat(glmfit.beta);
            
            %Analysis 2: fit on opto trials and fit on control trials and compare weights to see differences: (Hwang et al. 2018)
            trialidx_opto_test      = find(X_photoStim);
            temp                    = find(~X_photoStim);
            trialidx_ctrl_test      = temp(randperm(length(temp),length(trialidx_opto))); %take equal amount of trials in control condition
            
            trialidx_ctrl_train     = ~ismember(1:size(temptrialData.session_ID,1),[trialidx_opto_test; trialidx_ctrl_test]); %rest is training data
            
            %Actual fit without crossvalidation for control trials:
            glmfit                  = glmnet(allX(trialidx_ctrl_train,:),Y(trialidx_ctrl_train),'multinomial',glmoptions);
            
            Y_hat_choice            = glmnetPredict(glmfit,allX(trialidx_ctrl_test,:),glmoptions.lambda,'class');
            CV_ctrl_2(iMou)         = sum(Y_hat_choice == Y(trialidx_ctrl_test)) / numel(Y(trialidx_ctrl_test));
            
            Y_hat_choice            = glmnetPredict(glmfit,allX(trialidx_opto_test,:),glmoptions.lambda,'class');
            CV_opto_2(iMou)         = sum(Y_hat_choice == Y(trialidx_opto_test)) / numel(Y(trialidx_opto_test));
            
    end
end

% nModels = size(B,1); %count number of different models that were fit
% CV_Perf = CV_Perf(1:nModels,:); %trim the cross validated performance

%% 
Colors  = parula(nMice);
Colors = getPyPlot_cMap('brg',nMice);

%% Analysis 1: Make figure of performance of the model trained on control or opto trials only:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.35 .08 .27]); hold all;

for iMou = 1:nMice
%     plot([1 2],[CV_ctrl_1(iMou) CV_opto_1(iMou)],'.-','Color','k','MarkerSize',30); hold all;
    plot([1 2],[CV_ctrl_1(iMou) CV_opto_1(iMou)],'.-','Color',Colors(iMou,:),'MarkerSize',30); hold all;
end

% p = signrank(CV_ctrl_1,CV_opto_1);
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions

%Bayesian statistics:
tempbf      = bf.ttest(CV_ctrl_1,CV_opto_1);
bfsymb      = MOL_BFtoSymbol(tempbf);
tempd       = computeCohen_d(CV_ctrl_1,CV_opto_1,'paired');
fprintf('Bayesian ttest (cv correct trained control or trained opto): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',size(CV_ctrl_1,1),tempd,tempbf)
bayesstar([1 2],tempbf);

%Make up:
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0.33 0.75])
xlim([0.7 2.3])
ylabel('Cross-validated % Correct')
% legend(legendhandles,uMice,'Location','NorthWest','FontSize',8); legend boxoff;
% text(1,0.5,sprintf('p=%1.3f',p),'FontSize',15)

%% Analysis 1: Make figure of weights of the history and non-history weights when trained on control or opto trials only:
idx_history         = [3 4 5];
idx_nonhistory      = ~ismember(1:7,idx_history);

figure; set(gcf,'color','w','units','normalized','Position', [0.3 0.35 .08 .27]); hold all;

for iMou = 1:nMice
    temp_ctrl = B_ctrl_1{iMou}(idx_history,:);
    temp_opto = B_opto_1{iMou}(idx_history,:);
    B_sum_hist_ctrl(iMou) = sum(abs(temp_ctrl(:)));
    B_sum_hist_opto(iMou) = sum(abs(temp_opto(:)));
    plot([1 2],[B_sum_hist_ctrl(iMou) B_sum_hist_opto(iMou)],'.-','Color',Colors(iMou,:),'MarkerSize',30); hold all;
    
end

% p = signrank(B_sum_hist_ctrl,B_sum_hist_opto);
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions

%Bayesian statistics:
tempbf      = bf.ttest(B_sum_hist_ctrl,B_sum_hist_opto);
bfsymb      = MOL_BFtoSymbol(tempbf);
tempd       = computeCohen_d(B_sum_hist_ctrl,B_sum_hist_opto,'paired');
fprintf('Bayesian ttest (sum of weights history control vs opto): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',size(B_sum_hist_ctrl,1),tempd,tempbf)
bayesstar([1 2],tempbf);

%Make up:
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
% ylim([0.33 0.75])
xlim([0.7 2.3])
ylabel('a.u.')
% legend(legendhandles,uMice,'Location','NorthWest','FontSize',8); legend boxoff;
% text(1,4,sprintf('p=%1.3f',p),'FontSize',15)

%% Analysis 2: Make figure of performance of the model when crossvalidated on control or on opto trials:
figure; set(gcf,'color','w','units','normalized','Position', [0.5 0.35 .08 .27]); hold all;

for iMou = 1:nMice
    plot([1 2],[CV_ctrl_2(iMou) CV_opto_2(iMou)],'.-','Color',Colors(iMou,:),'MarkerSize',30); hold all;
end

% p = signrank(CV_ctrl_2,CV_opto_2);
% sigstar([1 2],p) %use sigstar function to identify signficantly different conditions

%Bayesian statistics:
tempbf      = bf.ttest(CV_ctrl_2,CV_opto_2);
bfsymb      = MOL_BFtoSymbol(tempbf);
tempd       = computeCohen_d(CV_ctrl_2,CV_opto_2,'paired');
fprintf('Bayesian ttest (cv correct trained control and validated on opto): %d sessions, Cohen''s d = %1.3f, BF10=%3.3f\n',size(CV_ctrl_2,1),tempd,tempbf)
bayesstar([1 2],tempbf);

%Make up:
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0.33 0.75])
xlim([0.7 2.3])
ylabel('Cross-validated % Correct')
% legend(legendhandles,uMice,'Location','NorthWest','FontSize',8); legend boxoff;
% text(1,0.5,sprintf('p=%1.3f',p),'FontSize',15)


%% 


%% Loop over mice and then sessions:

mouseids                = unique(sessionData.mousename)';
nMice                   = length(mouseids);
B                       = {};
CV_Perf                 = NaN(10,length(mouseids));
% Y_Hat                   = NaN(10,length(mouseids),5000);

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
    
    % Choice history model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Choice history';
    X                           = [X_LastLick];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    %Log Sensory model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory';
    X                           = [X_logvisualOriChange X_logaudioFreqChange];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % Sensory + History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory + Choice History';
    X                           = [X_trialNumber X_logvisualOriChange X_logaudioFreqChange X_LastLick];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % FULL: Sensory Choice History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Full';
    X                           = [X_trialNumber X_LastLick [0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)] [0; X_logvisualOriChange(1:end-1)] [0; X_logaudioFreqChange(1:end-1)] X_logvisualOriChange X_logaudioFreqChange];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % Sensory Choice History model:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Sensory and Choice history + opto';
    X                           = [X_trialNumber X_logvisualOriChange X_logaudioFreqChange X_LastLick X_photoStim];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    
    % FULL model with opto:
    iModel                      = iModel + 1;
    Model_string{iModel}        = 'Full+opto';
    X                           = [X_trialNumber X_LastLick [0; X_visualCorrect(1:end-1)] [0; X_audioCorrect(1:end-1)] [0; X_logvisualOriChange(1:end-1)] [0; X_logaudioFreqChange(1:end-1)] X_logvisualOriChange X_logaudioFreqChange X_photoStim];
    [B{iModel,iMou},CV_Perf(iModel,iMou),~] = MOL_fitGLM_session(X,Y_Response,params);
    

    %     Export data to .csv file
    %     X                   = [X_visualOriChange X_audioFreqChange X_visualCorrect X_audioCorrect];
    %     rootdir             = 'E:\Documents\PhD\TempExportPietro\';
    %     csvwrite(sprintf('%s%s.csv',rootdir,mouseids{iMou}),[Y_Response X],0) %writes matrix M to file filename as comma-separated values.
    
end

nModels = size(B,1); %count number of different models that were fit
CV_Perf = CV_Perf(1:nModels,:); %trim the cross validated performance

%% Analysis 4: Model #7 compared to #9:
figure; set(gcf,'color','w','units','normalized','Position', [0.7 0.35 .08 .27]); hold all;
ctrl = 7;
opto = 9;

for iMou = 1:nMice
%     plot([1 2],[CV_ctrl_1(iMou) CV_opto_1(iMou)],'.-','Color','k','MarkerSize',30); hold all;
    plot([1 2],[CV_Perf(ctrl,iMou) CV_Perf(opto,iMou)],'.-','Color',Colors(iMou,:),'MarkerSize',30); hold all;
end

%Bayesian statistics:
tempbf      = bf.ttest(CV_Perf(ctrl,:),CV_Perf(opto,:));
bfsymb      = MOL_BFtoSymbol(tempbf);
tempd       = computeCohen_d(CV_Perf(ctrl,:),CV_Perf(opto,:),'paired');
fprintf('Bayesian ttest (Model #%d vs #%d): Cohen''s d = %1.3f, BF10=%3.3f\n',ctrl,opto,tempd,tempbf)
bayesstar([1 2],tempbf);

%Make up:
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0.33 0.75])
xlim([0.7 2.3])
ylabel('Cross-validated % Correct')

%% Analysis 5: Model #8 compared to #10:
figure; set(gcf,'color','w','units','normalized','Position', [0.8 0.35 .08 .27]); hold all;
ctrl = 8;
opto = 10;

for iMou = 1:nMice
%     plot([1 2],[CV_ctrl_1(iMou) CV_opto_1(iMou)],'.-','Color','k','MarkerSize',30); hold all;
    plot([1 2],[CV_Perf(ctrl,iMou) CV_Perf(opto,iMou)],'.-','Color',Colors(iMou,:),'MarkerSize',30); hold all;
end

%Bayesian statistics:
tempbf      = bf.ttest(CV_Perf(ctrl,:),CV_Perf(opto,:));
bfsymb      = MOL_BFtoSymbol(tempbf);
tempd       = computeCohen_d(CV_Perf(ctrl,:),CV_Perf(opto,:),'paired');
fprintf('Bayesian ttest (Model #%d vs #%d): Cohen''s d = %1.3f, BF10=%3.3f\n',ctrl,opto,tempd,tempbf)
bayesstar([1 2],tempbf);

%Make up:
set(gca, 'XTick', 1:2,'XTickLabels', {'Control' 'Photostim'},'XTickLabelRotation',45)
ylim([0.33 0.75])
xlim([0.7 2.3])
ylabel('Cross-validated % Correct')

