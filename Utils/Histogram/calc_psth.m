function [hist_mat] = calc_psth(events_ts,spikes_ts,params)
% Standard function to compute firing rates over time aligned to events
% input:
% - events_ts: timestamps for specific events
% - timestamps of all spikes
% - parameters (Parameter settings are to be defined elsewhere)
% output:
% psth matrix of size (M by N) where M is number of events/trials, and N is
% number of timebins
% Matthijs oude Lohuis 2017-2021

%% Get spiking times:
trial_ts            = cell(1,length(events_ts));
for ev=1:length(events_ts) % Get spikes within window around event ev:
    trial_ts{ev}     = spikes_ts(spikes_ts>events_ts(ev)+params.t_pre & spikes_ts<=events_ts(ev)+params.t_post)-events_ts(ev);
end

%% Define time axis:
if ~isfield(params,'edges') || ~isfield(params,'xtime')
    params.edges                = params.t_pre:params.binsize:params.t_post;                    %Define edges of the bins
    params.xtime                = [params.t_pre:params.binsize:params.t_post-params.binsize] + params.binsize/2;     %Define time axis
    params.nTimebins            = length(params.xtime);
end

%% Histogram:
hist_mat = NaN(length(events_ts),params.nTimebins);
for ev=1:length(events_ts)
    temp                = histc(trial_ts{ev},params.edges) * 1e6/params.binsize;
    if ~isempty(temp) %No spikes whatsover in the trial
        hist_mat(ev,:)      = temp(1:end-1); %Store + discard last value (=exactly values at endbin value)
    end
end

%% Apply smoothing if requested:
if params.smoothing
    if strcmp(params.conv_win,'flat') %Construct smoothing window:
        win                 = ones(1,params.conv_twin/params.binsize)/(params.conv_twin/params.binsize);
    elseif strcmp(params.conv_win,'gaussian')
        N                   = params.conv_twin/params.binsize;
        alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
        win                 = gausswin(N,alpha); %convolution with gaussian
        win                 = win/sum(win); %normalized
    elseif strcmp(params.conv_win,'chg') %causal half-gaussian:
        N                   = params.conv_twin/params.binsize;
        alpha               = ((N-1)/(params.conv_sigma/params.binsize))/2; %? = (N – 1)/(2?)
        win                 = gausswin(N,alpha); %convolution with gaussian
        win(1:round(N/2)-1) = 0; %set acausal bins to zero;
%         win((round(N/2)-1):end) = 0; %set acausal bins to zero;
        win                 = win/sum(win); %normalized
    end
    %Smooth the individual trials:
    hist = cell(1,length(events_ts));
    for ev = 1:length(events_ts)
        hist{ev}                = padarray(hist_mat(ev,:),[0 round(length(win)/2)],'symmetric','both'); %pad the array on both sides for convolution
        hist{ev}                = conv(hist{ev},win,'valid'); %Take only the valid overlapping center of convolution
        hist_mat(ev,:)          = hist{ev}(1:params.nTimebins); %slight correction to get same size (edge vs convolution)
    end
end

%% Z-SCORE
if params.zscore
    %Get baseline period:
    bsl_idx = params.xtime>params.twin_baseline_start & params.xtime<=params.twin_baseline_stop;
    %Get std per trial baseline period:
    std_psth_bsl_win = nanstd(hist_mat(:,bsl_idx),0,2);
    std_psth_bsl_win = std_psth_bsl_win(~isnan(std_psth_bsl_win));
    %Each trial minus its mean baseline activity divided by mean standard
    %deviation during baseline all trials:
    if nanmean(std_psth_bsl_win)==0
        warning('zero spikes found during baseline, no z-score computed, z set to zero')
        hist_mat(:,:) = 0;
    else
        for ev = 1:size(hist_mat,1)
            hist_mat(ev,:) = (hist_mat(ev,:) - nanmean(hist_mat(ev,bsl_idx)))/nanmean(std_psth_bsl_win);
        end
    end
end

end
