function filteredspikeData = MOL_filterNeurons(varargin)
%% Get input arguments:
tempsessionData         = varargin{1};
temptrialData           = varargin{2};
tempspikeData           = varargin{3};

%%
fprintf('Filtering Neurons... \n')

%% Criteria:
criteria.max_ISI_FA     = 0.01; % Select on fraction of spikes within absolute refractory period (1.5ms)
criteria.min_ID         = 5;    % Select on Isolation Distance
criteria.min_avg_fr     = 0.2;  % Select on minimum average firing rate throughout the session
criteria.min_coverage   = 0.7;  % Select on coverage throughout the session (fraction of session):
criteria.sign_resp      = 0;    % Whether any condition needs to show a significant response

if criteria.sign_resp
    if strfind(tempsessionData.Experiment{1},'Bars')
        params                      = params_histresponse_pmal(); % All time is in microseconds
        params.trialrespmethod      = 'total';      %Whether to bin (and smooth) individual trials 'individual' or all trials together ('total')
        params.respcalcmethod       = 'max';        %Which method to calculate response {'max','mean','div','AUC'}
        params.zscore               = 1;
        params.minzscore            = 2;
        params.eventofinterest      = 'stimStart';
        params.trialcondidx = [];
        for iOri = unique(temptrialData.visualOri)'
            for iSpeed = unique(temptrialData.visualSpeed)'
                idx = temptrialData.visualOri == iOri & temptrialData.visualSpeed == iSpeed  & temptrialData.hasphotostim == 0;
                params.trialcondidx = [params.trialcondidx idx];
            end
        end
        params.trialcondidx = params.trialcondidx==1;
    else    params                  = params_histresponse(); % All time is in microseconds
        params.eventofinterest      = 'stimChange';
        params.condofinterest{1}    = 'trialType';
    end
    
end

%% Initialize selected index:
criteria.idx = true(size(tempspikeData.ch));

%% Select on Refractory period violations:
if criteria.max_ISI_FA && isfield(tempspikeData,'QM_ISI_FA')
    ISI_FA              = tempspikeData.QM_ISI_FA<criteria.max_ISI_FA;
    fprintf('Excluded %d/%d: Refractory period violation rate too high\n',sum(~ISI_FA & criteria.idx),sum(criteria.idx))
    criteria.idx    = criteria.idx & ISI_FA;
end

%% Select on Isolation Distance:
if criteria.min_ID && isfield(tempspikeData,'QM_IsolationDistance')
    ID              = tempspikeData.QM_IsolationDistance>criteria.min_ID;
    fprintf('Excluded %d/%d: too low Isolation Distance\n',sum(~ID & criteria.idx),sum(criteria.idx))
    criteria.idx    = criteria.idx & ID;
end

%% Select on minimum average firing rate throughout the session:
if criteria.min_avg_fr && isfield(tempspikeData,'avg_fr')
    AVGFR = tempspikeData.avg_fr>criteria.min_avg_fr;
    fprintf('Excluded %d/%d: too low average firing rate throughout session\n',sum(~AVGFR & criteria.idx),sum(criteria.idx))
    criteria.idx = criteria.idx & AVGFR;
end

%% Select on coverage throughout the session (fraction of session):
if criteria.min_coverage && isfield(tempspikeData,'coverage')
    COV = tempspikeData.coverage>criteria.min_coverage;
    fprintf('Excluded %d/%d: too little coverage of the session\n',sum(~COV & criteria.idx),sum(criteria.idx))
    criteria.idx = criteria.idx & COV;
end

%% Select only those cells that have a significant response to one condition:
if criteria.sign_resp
    for iNeuron = 1:length(tempspikeData.ch)
        sign_resp_idx(iNeuron,1) = 0;
        %Get histogram per neuron:
        events_ts = temptrialData.(params.eventofinterest)(temptrialData.session_ID == tempspikeData.session_ID(iNeuron));
        spikes_ts = tempspikeData.ts{iNeuron};
        
        [edges,hist_mat]    = calc_psth(events_ts,spikes_ts,params);
        
        for iCond = 1:size(params.trialcondidx,2)
            if strcmp(tempsessionData.Experiment{1},'Bars') %Get correct time windows:
                params.twin_resp_stop_trial = temptrialData.stimEnd(params.trialcondidx(:,iCond)==1) - temptrialData.stimStart(params.trialcondidx(:,iCond)==1);  %different response window per trial (variable stimulus length)
            end
            %Get responses:
            [resp,~]                        = calc_resp_from_psth(edges,hist_mat(params.trialcondidx(:,iCond),:),params);
            if mean(resp)>params.minzscore
                sign_resp_idx(iNeuron,1) = 1;
            end
        end
    end
    SIGN = sign_resp_idx;
    fprintf('Excluded %d/%d: no significant response to any condition\n',sum(~SIGN & criteria.idx),sum(criteria.idx))
    criteria.idx = criteria.idx & SIGN;
end

%% Filter neurons based on resulting index:

datafields = fieldnames(tempspikeData);
filteredspikeData = struct();
for field = 1:length(datafields)
    filteredspikeData.(datafields{field}) = tempspikeData.(datafields{field})(criteria.idx);
end

fprintf('Subselected %d neurons out of %d\n',length(filteredspikeData.ch),length(tempspikeData.ch))

end