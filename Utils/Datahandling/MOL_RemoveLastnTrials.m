function trialData = MOL_RemoveLastnTrials(trialData,lastn)

uniqueses = unique(trialData.session_ID);
for iSes = uniqueses'
    trialfieldnames = fieldnames(trialData);
    nTotalTrials = max(trialData.trialNum(strcmp(trialData.session_ID,iSes)));
    selec = trialData.trialNum > nTotalTrials - lastn & strcmp(trialData.session_ID,iSes);
    if sum(selec) ~= lastn; error('more or less trials to be discarded than requested'); end
    
    for iField = 1:length(trialfieldnames)
        trialData.(trialfieldnames{iField}) = trialData.(trialfieldnames{iField})(~selec);
    end
end

fprintf('Trimmed all sessions: last %2.0f trials removed\n\n',lastn)
end