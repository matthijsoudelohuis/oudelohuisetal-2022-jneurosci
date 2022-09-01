%% Script that analyzes primary behavioral measures of performance across the three task versions:
% MOL (C) 2021

%% 
params.Experiments          = {'ChangeDetectionConflictDecor' {'BehaviorConflict' 'BehaviorPsychophysics'}};
params.ExperimentLabels     = {'NE' 'MST'};
params.nExperiments         = length(params.Experiments);

params                      = MOL_getColors_CHDET(params);

%% General figure settings:
set(0,'defaultAxesFontSize',20)
set(0,'Defaultlinelinewidth',5)
set(0,'DefaultAxesFontName','Arial')

%% Get the data:
sessionData = struct(); trialData = struct();
for iExp = 1:length(params.Experiments)
    [Data]                  = MOL_GetData('E:','CHDET',params.Experiments{iExp},{},{},{'sessionData' 'trialData'});
    sessionData             = AppendStruct(sessionData,Data.sessionData);
    trialData               = AppendStruct(trialData,Data.trialData);
end
trialData               = MOL_RemoveLastnTrials(trialData,20); %Remove last 20 trials

%% Fit each session:
nSessions               = length(sessionData.mousename);

FullParam               = NaN(nSessions,8); %Init output variable

for iSes = 1:nSessions %Fit each session
    fprintf('\nFitting session %d/%d \n\n',iSes,nSessions);
    %Get the data for this session only:
    [tempsessionData,temptrialData]                 = MOL_getTempPerSes(sessionData.session_ID(iSes),sessionData,trialData);
    %Get response rates per condition:
    [visconditions,auconditions,FullRespMat,~]      = MOL_Psy_GetTrialMat(tempsessionData,temptrialData);
    %Correct dimensions for some sessions:
    if numel(visconditions)==5 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end,:,:)];
        auconditions = [2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==6 && numel(auconditions)==4
        FullRespMat = [FullRespMat; FullRespMat(end-1:end,:,:)];
        auconditions = [1/256 2/256 8/256 32/256 64/256 128/256];
    end
    if numel(visconditions)==4 && numel(auconditions)==3
        FullRespMat = [FullRespMat(1:2,:,:); FullRespMat(2:end,:,:)];
        auconditions = [2/256 8/256 32/256 128/256];
    end
    if numel(visconditions)==5 && numel(auconditions)==8
        idx = [1 3 5 7 8];
        FullRespMat = FullRespMat([1 idx+1],:,:);
        auconditions = auconditions(idx);
    end
    
    %Rewrite fullresponse table to contingency table for model fitting (without conflict trials)
    ctable              = NaN(3,3,numel(visconditions));
    ctable(1,:,:)       = permute(FullRespMat(2:end,1,:),[3 1 2]);              %Auditory
    ctable(2,:,:)       = permute(FullRespMat(1,2:end,:),[1 3 2]);              %Visual
    ctable(3,:,:)       = repmat(permute(FullRespMat(1,1,:),[1 3 2]),1,1,numel(visconditions));    %Probe trials
    
    %Align visual and auditory intensities
    %visual and auditory conditions should be normalized such that
    %value of 1 corresponds to expected asymptotic dmax
    normconditions = [auconditions / max(auconditions); visconditions / max(visconditions)];
    
    %Fit psychometric 2ADC model:
    %     fprintf('Fitting session %d/%d, of animal %d/%d\n',ses,length(sesselec),mou,length(mouseids));
    [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc_editMOL(ctable,normconditions,[]);
    
    %Store parameters:
    %theta_est = 8 parameters: % 3 d' -related parameters and one 'c' parameter for each location
    % order is daud, dvis, naud, nvis, s50aud, s50vis, caud, vis
    % dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M); c = theta(3*M+1:end);
    theta_est(5) = theta_est(5) * max(auconditions); %undo normalization of conditions
    theta_est(6) = theta_est(6) * max(visconditions); %undo normalization of conditions
    FullParam(iSes,:)            = theta_est; %store parameters
    
    %     theta_err(5) = theta_err(5) * max(auconditions); %undo normalization of conditions
    %     theta_err(6) = theta_err(6) * max(visconditions); %undo normalization of conditions
    %     FullParamErr(iSes,:)            = theta_err; %store parameters
    
end
% 
% %% Plot average rates for each cohort:
% for iExp = 1:params.nExperiments
%     
%     figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.1 .7 .7]); hold all;
%     
%     expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
%     
%     for iAnimal = 1:length(expanimals)
%         
%         [tempsessionData,~]        = MOL_getTempPerSes(sessionData.session_ID(strcmp(sessionData.mousename,expanimals{iAnimal})),sessionData,trialData);
%         
%         idx_ses  = strcmp(sessionData.mousename,expanimals{iAnimal});
%         meantheta_est             = median(FullParam(idx_ses,:),1);
%         
%         % Settings:
%         switch tempsessionData.auChangeUnit{1}
%             case 'Hz'
%                 params.auprobepos = 0.5;
%                 params.auticks = [100 4000];
%                 params.auticklabels = ['Probe' num2cell(params.auticks)];
%                 params.auxaxislabel  = 'Delta frequency (Hz)';
%                 params.auystatslabel = 'Auditory threshold (Hz)';
%             case 'Oct'
%                 params.auprobepos = 0.001;
%                 params.auticks = [1/256 1/64 1/8 1/2];
%                 params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
%                 params.auxaxislabel  = 'Delta frequency (Oct)';
%                 params.auystatslabel = 'Auditory threshold (partial octave)';
%         end
%         
%         params.visprobepos     = 0.5;
%         params.visticks        = [5 30 90];
%         params.vistickslabels  = ['Probe' num2cell(params.visticks)];
%         params.visxaxislabel   = 'Delta orientation (Degrees)';
%         params.visystatslabel  = 'Visual threshold (Degrees)';
%         
%         params.yticks          = [0 0.25 0.5 0.75 1];
%         
%         % Generate contingency table from fitted parameters:
%         [~,xvals_fit_vis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(meantheta_est,params);
%         %Make common x-axis for the auditory to merge from both cohorts:
%         xvals_fit_au = 10.^(linspace(log10(0.001),log10(1/2),1000)); %logarithmic spacing
%         
%         %Audio:
%         subplot(1,2,1); hold all;
%         plot(xvals_fit_au,squeeze(ctable_fit_mat(1,1,:)),'-','Color',[1 0 0],'LineWidth',3);
%         plot(xvals_fit_au,squeeze(ctable_fit_mat(1,2,:)),':','Color',[0.2 0.2 1],'LineWidth',3);
%         
%         %Visual:
%         subplot(1,2,2); hold all;
%         plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,2,:)),'-','Color',[0 0 1],'LineWidth',3);
%         plot(xvals_fit_vis,squeeze(ctable_fit_mat(2,1,:)),':','Color',[1 0.2 0.2],'LineWidth',3);
%     end
%     
%     params.auprobepos = 0.001;
%     params.auticks = [1/256 1/64 1/8 1/2];
%     params.auticklabels = {'Probe' '1/256' '1/64' '1/8' '1/2'};
%     params.auxaxislabel  = 'Delta frequency (Oct)';
%     params.auystatslabel = 'Auditory threshold (partial octave)';
%     
%     MOL_Psy2Sided_FigMakeup(params)
%     
% end

%% Show Dprime for all experiments
%Get all the dprimes in one multidimensional variable:
datatoplot      = NaN(params.nExperiments,2,20,20); %init output var (dim = experiments, modality, mouse, session)
for iExp = 1:params.nExperiments
    expanimals                 = unique(sessionData.mousename(ismember(sessionData.Experiment,params.Experiments{iExp})));
    
    for iAnimal = 1:length(expanimals)
        sesidx                  = ismember(sessionData.mousename,expanimals(iAnimal));
        datatoplot(iExp,1,iAnimal,1:sum(sesidx))        = FullParam(sesidx,1);
        datatoplot(iExp,2,iAnimal,1:sum(sesidx))        = FullParam(sesidx,2);
    end
end

datatoplot_re       = reshape(datatoplot,params.nExperiments,2,20*20); %reshape
% mediandatatoplot    = nanmedian(datatoplot,4); %average over sessions
meandatatoplot      = nanmean(datatoplot_re,3); %average over sessions
errortoplot         = nanstd(datatoplot_re,[],3); %average over sessions
nSessions           = sum(~isnan(datatoplot_re),3); %average over sessions
errortoplot = errortoplot ./ sqrt(nSessions);
errortoplot = errortoplot ./ sqrt(nSessions);

figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .2 .4]); hold all;
% params.offset       = 1; %Offset of bars relative to each other

for iExp= 1:params.nExperiments
    for iMod = 1:2
%         errorbar((iExp-1)*2+iMod,meandatatoplot(iExp,iMod),errortoplot(iExp,iMod),'.','MarkerSize',20,'LineWidth',0.01,'Color','k')
        %With std across all sessions:
        errorbar((iExp-1)*2+iMod,meandatatoplot(iExp,iMod),errortoplot(iExp,iMod),'.','MarkerSize',0.1,'LineWidth',15,'Color',params.colors_modalities{iMod})
        errorbar((iExp-1)*2+iMod,meandatatoplot(iExp,iMod),errortoplot(iExp,iMod),'.','MarkerSize',25,'LineWidth',0.01,'Color','k')
%         errorbar(1+params.offset*(iExp-1)+params.offset/3*(iMod-1),nanmean(meandatatoplot(iExp,iMod,:),3),nanstd(meandatatoplot(iExp,iMod,:),[],3),'k.','MarkerSize',1,'LineWidth',5)
    end
end

ylabel('Dprime')
set(gca,'XTick',1:4,'XTickLabels',{'Aud' 'Vis' 'Aud' 'Vis'},'XTickLabelRotation',45)
ylim([-0.2 3])
xlim([0.5 4.5])

%Statistics: (Comparing audio d-prime for unisensory vs multisensory)
Xdata = datatoplot_re(1,1,:); Xdata = Xdata(:);
Ydata = datatoplot_re(2,1,:); Ydata = Ydata(:);
tempBF = bf.ttest(Xdata,Ydata);
tempCohD = computeCohen_d(Xdata,Ydata,'independent');
fprintf('Naive vs trained auditory dprime: Bayesian t-test (n=%d naive sessions, n=%d, trained sessions): Cohen''s d: %1.2f BF=%1.2e\n',sum(~isnan(Xdata)),sum(~isnan(Ydata)),tempCohD,tempBF)

bayesstar([1 3],tempBF)

%Statistics: (Comparing visual d-prime for unisensory vs multisensory)
Xdata = datatoplot_re(1,2,:); Xdata = Xdata(:);
Ydata = datatoplot_re(2,2,:); Ydata = Ydata(:);
tempBF = bf.ttest(Xdata,Ydata);
tempCohD = computeCohen_d(Xdata,Ydata,'independent');
fprintf('Naive vs trained visual dprime: Bayesian t-test (n=%d naive sessions, n=%d, trained sessions): Cohen''s d: %1.2f BF=%1.2e\n',sum(~isnan(Xdata)),sum(~isnan(Ydata)),tempCohD,tempBF)

bayesstar([2 4],tempBF)

