function [d1,d2,c1,c2] = MOL_2ADC_Full(varargin)

%% Get input arguments:
sessionData     = varargin{1};
trialData       = varargin{2};

%% General settings:
showIndFig          = 0;
showResFig          = 1;
set(0,'defaultAxesFontSize',20)
sessioncounter      = 0;

%% Remove all trials that are not maximal change (or probe):
if strcmp(sessionData.auChangeUnit(1),'Hz')
    idx         = (abs(trialData.visualOriChange)==90 | isnan(trialData.visualOriChange)) & (abs(trialData.audioFreqChange)==4000 | isnan(trialData.audioFreqChange));
else strcmp(sessionData.auChangeUnit(1),'Oct')
    idx         = (abs(trialData.visualOriChange)==90 | isnan(trialData.visualOriChange)) & (abs(trialData.audioOctChange)==0.5 | isnan(trialData.audioFreqChange));
end
datafields = fieldnames(trialData);
for field = 1:length(datafields)
    trialData.(datafields{field}) = trialData.(datafields{field})(idx);
end

%% Loop over mice and then sessions:
mouseids = unique(sessionData.mousename)';
for mou = 1:length(mouseids)
    %Get the relevant sessions for this mouse:
    sesselec        = unique(sessionData.session_ID(strcmp(sessionData.mousename,mouseids{mou})))';
    
    for ses = 1:length(sesselec) %Loop over sessions for this mouse
        sessioncounter = sessioncounter+1;         %Add one to the overall counter
        sesid = sesselec(ses);                      %Get sessionID for this session to subselect trials and info about this sessions:
        [tempsessionData,temptrialData] = MOL_getTempPerSes(sesid,sessionData,trialData);%Get the sessionData for each session individually
        [dVis(mou,ses),dAud(mou,ses),cVis(mou,ses),cAud(mou,ses)] = MOL_Fit_2ADC_Full_Session(tempsessionData,temptrialData,showIndFig); %Compute 2ADC dprime
    end
    
end

%% 
dVis(dVis==0) = NaN;
dAud(dAud==0) = NaN;

%%

if showResFig
    f = figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.4 .8 .4]);
    subplot(1,2,1);
    plot(1:length(mouseids),dVis,'b.','MarkerSize',30); hold all;
    
    mVis = nanmean(dVis(:));
    stdVis = nanstd(dVis(:)) ;%/sqrt(sum(~isnan(dVis(:))));
    plot(length(mouseids)+1,mVis,'k.','MarkerSize',30)
    errorbar(length(mouseids)+1,nanmean(dVis(:)),stdVis,stdVis,'k','LineWidth',2)
    
    ylabel('d-Prime Visual')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 4])
    xlim([0.5 length(mouseids)+1.5])
    
    subplot(1,2,2);
    plot(1:length(mouseids),dAud,'r.','MarkerSize',30); hold all;
    
    mAud = nanmean(dAud(:));
    stdAud = nanstd(dAud(:)); %/sqrt(sum(~isnan(dAud(:))));
    plot(length(mouseids)+1,mAud,'k.','MarkerSize',30)
    errorbar(length(mouseids)+1,nanmean(dAud(:)),stdAud,stdAud,'k','LineWidth',2)
    
    ylabel('d-Prime Auditory')
    set(gca, 'XTick', 1:length(mouseids),'XTickLabels', mouseids,'XTickLabelRotation',45)
    ylim([0 4])
    xlim([0.5 length(mouseids)+1.5])
    
end
