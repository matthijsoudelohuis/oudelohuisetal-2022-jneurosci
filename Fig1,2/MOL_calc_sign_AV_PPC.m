function [sessionData,trialData,spikeData] = MOL_calc_sign_AV_PPC(sessionData,trialData,spikeData)
%% General parameters:
params                      = params_histresponse(); % All time is in microseconds
params.zscore               = 0;
params.smoothing            = 0;

params.AlignOn                = 'stimChange';      %On which timestamp to align as t=0

% parameters for window of interest:
params.t_start_0              = -0.5e6;
params.t_stop_0               = -0.1e6;
params.t_start_1              = 0e6;
params.t_stop_1               = 0.5e6;

% params.t_window_dur         = params.t_post-params.t_pre;
params.ttestalpha             = 0.01;

%% Init output fields:
spikeData.sign_incr_aud      = NaN(size(spikeData.session_ID));
spikeData.sign_decr_aud      = NaN(size(spikeData.session_ID));
spikeData.sign_incr_vis      = NaN(size(spikeData.session_ID));
spikeData.sign_decr_vis      = NaN(size(spikeData.session_ID));

%%
params.nSplits = 2;

nNeurons                = length(spikeData.ts);
lastsesid               = []; %memory var to keep track of whether neuron comes from different session and design matrix needs to be reconstructed
fprintf('Computing significant response for neuron        \n');

for iNeuron = 1:nNeurons %Loop over all neurons:
    fprintf(repmat('\b', 1, numel([num2str(iNeuron-1) num2str(nNeurons)])+2));
    fprintf('%d/%d\n',iNeuron,nNeurons);
    
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
    
    baseline      = nanmean(hist_mat(:,params.xtime>params.t_start_0 & params.xtime<params.t_stop_0),2);
    response      = nanmean(hist_mat(:,params.xtime>params.t_start_1 & params.xtime<params.t_stop_1),2);
    
    spikeData.sign_incr_vis(iNeuron) = bf.ttest(baseline(splits{1}),response(splits{1}),'tail','left'); %Ttest response greater than baseline
    spikeData.sign_decr_vis(iNeuron) = bf.ttest(baseline(splits{1}),response(splits{1}),'tail','right'); %Ttest baseline greater than response
    spikeData.sign_incr_aud(iNeuron) = bf.ttest(baseline(splits{2}),response(splits{2}),'tail','left');
    spikeData.sign_decr_aud(iNeuron) = bf.ttest(baseline(splits{2}),response(splits{2}),'tail','right');
    
%     spikeData.sign_incr_vis(iNeuron) = ttest(baseline(splits{1}),response(splits{1}),'alpha',params.ttestalpha,'tail','left'); %Ttest response greater than baseline
%     spikeData.sign_decr_vis(iNeuron) = ttest(baseline(splits{1}),response(splits{1}),'alpha',params.ttestalpha,'tail','right'); %Ttest baseline greater than response
%     spikeData.sign_incr_aud(iNeuron) = ttest(baseline(splits{2}),response(splits{2}),'alpha',params.ttestalpha,'tail','left');
%     spikeData.sign_decr_aud(iNeuron) = ttest(baseline(splits{2}),response(splits{2}),'alpha',params.ttestalpha,'tail','right');
    
end

%     
%     
% %% For each session for each neuron compute sign resp to idx:
% for sesid = unique(sessionData.session_ID)'
%     %% Get the relevant data for each session individually:
%     [temptrialData,tempspikeData] = MOL_getTempPerSes(sesid,trialData,spikeData);
%     
%     %% Construct tensor:
%     events_ts       = temptrialData.(params.AlignOn);
%     nNeurons        = length(tempspikeData.ts);
%     nTimebins       = length(params.t_pre:params.binsize:params.t_post-params.binsize);
%     nTrials         = length(events_ts);
%     tensor          = NaN(nNeurons,nTimebins,nTrials);
%     for iNeuron = 1:nNeurons
%         %Get histogram per neuron
%         [edges,hist_mat]        = calc_psth(events_ts,tempspikeData.ts{iNeuron},params);
%         if params.smoothing 
%             tensor(iNeuron,:,:)     = hist_mat';
%         else
%             tensor(iNeuron,:,:)     = hist_mat' / (1e6/params.binsize); %Keep actual spike counts
%         end
%     end
%     
%     %% Compute response and baseline during specified timewindow:
%     baseline    = NaN(nNeurons,nTrials);
%     resp_1      = NaN(nNeurons,nTrials);
%     resp_2      = NaN(nNeurons,nTrials);
%     
%     for iNeuron = 1:nNeurons
%         baseline(iNeuron,:)       = squeeze(sum(tensor(iNeuron,edges>params.t_start_0 & edges<params.t_stop_0,:)));
%         resp_1(iNeuron,:)         = squeeze(sum(tensor(iNeuron,edges>params.t_start_1 & edges<params.t_stop_1,:)));
%         resp_2(iNeuron,:)         = squeeze(sum(tensor(iNeuron,edges>params.t_start_2 & edges<params.t_stop_2,:)));
%     end
%     
%     %% Compute significant response:
%     sestrialidx                 = idx(strcmp(trialData.session_ID,sesid));
%     sign_firstbump              = NaN(nNeurons,1);
%     sign_secondbump              = NaN(nNeurons,1);
%     
%     for iNeuron = 1:nNeurons
%         sign_firstbump(iNeuron) = ttest(baseline(iNeuron,sestrialidx),resp_1(iNeuron,sestrialidx),'alpha',params.ttestalpha);
%         sign_secondbump(iNeuron) = ttest(baseline(iNeuron,sestrialidx),resp_2(iNeuron,sestrialidx),'alpha',params.ttestalpha);
%     end
%     
%     spikeData.sign_firstbump(strcmp(spikeData.session_ID,sesid)) = sign_firstbump;
%     spikeData.sign_secondbump(strcmp(spikeData.session_ID,sesid)) = sign_secondbump;
%     
% end


%Convert to logical:
spikeData.sign_incr_vis = spikeData.sign_incr_vis>3;
spikeData.sign_decr_vis = spikeData.sign_decr_vis>3;
spikeData.sign_incr_aud = spikeData.sign_incr_aud>3;
spikeData.sign_decr_aud = spikeData.sign_decr_aud>3;

%Convert to logical:
% spikeData.sign_incr_vis = spikeData.sign_incr_vis==1;
% spikeData.sign_decr_vis = spikeData.sign_decr_vis==1;
% spikeData.sign_incr_aud = spikeData.sign_incr_aud==1;
% spikeData.sign_decr_aud = spikeData.sign_decr_aud==1;

fprintf('Calculated significant auditory or visual response for %d neurons\n\n', length(spikeData.sign_incr_vis))

end