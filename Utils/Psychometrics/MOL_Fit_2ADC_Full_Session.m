function [d1,d2,c1,c2] = MOL_Fit_2ADC_Full_Session(sessionData,trialData,showFig,varargin)
% MOL 2017
% In the 2-ADC model, the probability of response at
% each location (Y = i) for each stimulus event (X) can
% be derived from the structural model (Equation 1)
% and decision rule (Equation 2). We illustrate the case
% for p(Y = 1|X). The other cases may be similarly
% derived.
% The probability of response at location 1 is the
% combined probability that the decision variable at
% location 1 exceeds the choice criterion at that location
% and that its magnitude (over its choice criterion) is the
% larger of the two locations.
%
% Input variable is trialData. From this struct the variable
% outcome is the 3x3 stimulus response outcome matrix:
%            |              Response                |
%   Stimulus |  "1"       |  "2"       |  "Catch"')
%   ------------------------------------------------|
%   "1"      |            |            |            |
%   ---------+--------------------------------------|
%   "2"      |            |            |            |
%   ------------------------------------------------|
%   Absent   |            |            |            |
%   ------------------------------------------------|

% Output is the sensitivity (d-prime) for Stimulus 1, 2
% and criterion for Stimulus 1 and 2.

% Adapted from:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Sridharan Devarajan, Aug 2014
% dsridhar@stanford.edu; sridharan.d@gmail.com
%
% Citation:
%   Sridharan, D., Steinmetz, N.A., Moore, T. and Knudsen, E.I. (2014).
%       Distinguishing bias from sensitivity effects in multialternative
%       detection tasks. J. Vision 14 (9):16, 1-32. doi: 10.1167/14.9.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%


%% Compute stim-response contingencies from trialData

%Compute over all trials:
if ~isempty(varargin)
    outcome = varargin{1};
else
    if sessionData.VisualLeftCorrectSide
        HitVisual           = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 1 & trialData.noResponse ~= 1);
        ErrorVisual         = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 0 & trialData.noResponse ~= 1);
        FAVisual            = nansum(trialData.leftCorrect~=1 & trialData.rightCorrect~=1    & trialData.correctResponse == 0 & strcmp(trialData.responseSide,'L'));
        ErrorAudio          = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 0 & trialData.noResponse ~= 1);
        HitAudio            = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 1 & trialData.noResponse ~= 1);
        FAAudio             = nansum(trialData.leftCorrect~=1 & trialData.rightCorrect~=1    & trialData.correctResponse == 0 & strcmp(trialData.responseSide,'R'));
        MissVisual          = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 0 & trialData.noResponse == 1);
        MissAudio           = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 0 & trialData.noResponse == 1);
        CorrectRejection    = nansum(trialData.leftCorrect~=1  & trialData.rightCorrect~=1   & trialData.correctResponse == 1 & trialData.noResponse == 1);
    else
        %Compute over all trials:
        HitAudio            = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 1 & trialData.noResponse ~= 1);
        ErrorAudio          = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 0 & trialData.noResponse ~= 1);
        FAAudio             = nansum(trialData.leftCorrect~=1 & trialData.rightCorrect~=1    & trialData.correctResponse == 0 & strcmp(trialData.responseSide,'L'));
        ErrorVisual         = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 0 & trialData.noResponse ~= 1);
        HitVisual           = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 1 & trialData.noResponse ~= 1);
        FAVisual            = nansum(trialData.leftCorrect~=1 & trialData.rightCorrect~=1    & trialData.correctResponse == 0 & strcmp(trialData.responseSide,'R'));
        MissAudio           = nansum(trialData.leftCorrect==1                                & trialData.correctResponse == 0 & trialData.noResponse == 1);
        MissVisual          = nansum(trialData.rightCorrect==1                               & trialData.correctResponse == 0 & trialData.noResponse == 1);
        CorrectRejection    = nansum(trialData.leftCorrect~=1  & trialData.rightCorrect~=1   & trialData.correctResponse == 1 & trialData.noResponse == 1);
    end
    
    % Store in matrix:
    outcome  = [    HitVisual     ErrorAudio  MissVisual;
        ErrorVisual   HitAudio    MissAudio;
        FAVisual      FAAudio     CorrectRejection ];
end

% Clip outcomes by at least setting every category to 1, otherwise dprime values reach impossible values:
outcome(outcome==0) = 1;

%% Fit 2 ADC model:

datamat = zeros(5,5);
datamat(3:5,3:5) = outcome;

filtmat = zeros(5,5);
filtmat(3:5,3:5) = 1;

%% perform ML estimation
disp('Optimization starts...')
[theta_est, theta_err, data_fit] = MLE_4ADC(datamat, filtmat);

%% post-processing

% compute log-likelihood
nobs     = sum(datamat(:));
prob_fit = data_fit/nobs;
LLFm     = datamat.*log(prob_fit);
LLF      = sum(LLFm(:));

%%% display
disp(' ');
disp('---------------------')
disp('Observed data')
disp(datamat(3:5,3:5))

disp('Fitted data')
disp(floor(data_fit(3:5,3:5)))

disp(' ')
disp('ML estimated parameters')
fprintf('d1:%2.2f +/- %2.3f, d2:%2.2f +/- %2.3f \n', theta_est(3), theta_err(3), theta_est(4), theta_err(4));
fprintf('c1:%2.2f +/- %2.3f, c2:%2.2f +/- %2.3f \n'  , theta_est(7), theta_err(7), theta_est(8), theta_err(8));
disp('---------------------')

%% Sensitivity and bias parameters

% theta_est = round(theta_est,2);
% theta_est = double(int32(theta_est*100))/100;

d1 = theta_est(3);
d2 = theta_est(4);
c1 = theta_est(7);
c2 = theta_est(8);

%% Show the response counts for the actual and fitted data:
if showFig
    f = figure;
    set(gcf,'units','normalized','Position',[0.2 0.1 0.5 0.7],'color','w')
    
    % Create the uitable
    t = uitable(f,'Data',[datamat(3:5,3:5); floor(data_fit(3:5,3:5))],...
        'ColumnName',{'Resp 1','Resp 2','No Resp'},...
        'RowName',{'Stim 1','Stim 2','No Stim', 'FitStim 1','FitStim 2','Fit No Stim'},...
        'ColumnWidth',{60}, 'FontS',18);
    
    subplot(2,2,1) %,plot(3)
    pos = get(subplot(2,2,1),'position');
    title('table')
    xlabel('xlabel')
    ylabel('ylabel')
    set(subplot(2,2,1),'yTick',[])
    set(subplot(2,2,1),'xTick',[])
    
    set(t,'units','normalized')
    set(t,'position',pos)
    
    %%  Show fit in triple heat plot:
    
    % f = figure;
    % set(gcf,'units','normalized','Position',[0.1 0.4 0.8 0.4],'color','w')
    xaxisstart = -2; xaxisres =0.01; xaxisstop = 5;
    xaxis = xaxisstart:xaxisres:xaxisstop;
    
    c1axispos = (c1-xaxisstart) / xaxisres;
    c2axispos = numel(xaxis) - (c2-xaxisstart) / xaxisres;
    
    X1stim = [0 1 0];
    X2stim = [0 0 1];
    
    for signal = 1:3
        subplot(2,2,signal+1)
        %Set signal for the three possible conditions:
        X1 = X1stim(signal);
        X2 = X2stim(signal);
        
        mu      = [d1*X1 d2*X2];
        [x1,x2] = meshgrid(xaxis,xaxis);
        F       = mvnpdf([x1(:) x2(:)],mu);
        F       = reshape(F,length(xaxis),length(xaxis));
        colormap('hot')
        
        imagesc(flipud(F)); hold on;
        set(gca,'YTick',1:100:size(F,2),'YTickLabel',fliplr(xaxis(1:100:size(F,2))))
        set(gca,'XTick',1:100:size(F,2),'XTickLabel',xaxis(1:100:size(F,2)))
        xlabel('Stim Signal 1'); ylabel('Stim Signal 2');
        title(sprintf('Bivariate PDF Signal=%d',signal-1));
        
        plot([c1axispos c1axispos],[numel(xaxis) c2axispos],'LineWidth',2,'color','w'); hold on;
        plot([1 c1axispos],[c2axispos c2axispos],'LineWidth',2,'color','w'); hold on;
        plot([c1axispos numel(xaxis)],[c2axispos 1],'LineWidth',2,'color','w'); hold on;
        
    end
end



end