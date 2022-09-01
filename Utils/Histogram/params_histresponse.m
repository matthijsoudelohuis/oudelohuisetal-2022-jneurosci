function params = params_histresponse(varargin) 
% Parameter settings calculating responses and showing histogram
% All time is in microseconds

if nargin==0
    params                      = struct();
else params = varargin{1};
end

%% parameters for window of interest:
params.t_pre                = -1e6;
params.t_post               = 1.5e6;

%% parameters for binning:
%E.g. method can be either large squared bins or 1 ms bins with smoothing window ('smoothing') or combination thereof
params.binsize              = 1e3;          %Size of the bins
params.smoothing            = 1;    
params.conv_win             = 'gaussian';   %Type of window used for smoothing {flat, gaussian, chg)
params.conv_sigma           = 0.025e6;        %sd of gaussian window for smoothing
params.conv_twin            = round((9*params.conv_sigma)/params.binsize)*params.binsize;         %Window size for smoothing

%% Construct bin edges and time axis
params.edges                = [params.t_pre:params.binsize:params.t_post] - params.binsize/2;                    %#ok<NBRAK> Define edges of the bins
params.xtime                = params.t_pre:params.binsize:params.t_post-params.binsize;                        %Define time axis
params.nTimebins            = length(params.xtime); %number of time bins

%% Z-score?
params.zscore               = 1;            %Whether to convert psth matrix to zscore (x-meanxbaseline/sdbaseline)

%% parameters for calculating response
params.respcalcmethod       = 'mean';        %Which method to calculate response {'max','mean','div','AUC'}
params.trialrespmethod      = 'individual';
params.twin_baseline_start  = -1e6;         %Start window for calculating baseline
params.twin_baseline_stop   = -0.2e6;       %End window baseline
params.twin_resp_start      = 0e6;          %Start window for calculating response
params.twin_resp_stop       = 1e6;          %End window response %note that if stimulus is of variable length give input as 'twin_resp_stop_trial'
params.subtr_baseline       = 0;            %Subtract baseline or not

%% parameters for PSTH plot
params.plottype             = 'errorline';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 
% params.plottype             = 'bins';  %Either plot 'errorline' with mean and SEM, or plot 'bins' to plot classical bins 

end