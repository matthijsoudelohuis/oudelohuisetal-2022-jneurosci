%% example script for fitting an m-ADC (here, 4-ADC) model to data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Sridharan Devarajan, Aug 2014
% dsridhar@stanford.edu; sridharan.d@gmail.com
%
% Citation: 
%   Sridharan, D., Steinmetz, N.A., Moore, T. and Knudsen, E.I. (2014).
%       Distinguishing bias from sensitivity effects in multialternative 
%       detection tasks. J. Vision 14 (9):16, 1-32. doi: 10.1167/14.9.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% stimulus-response contingency matrix of response counts
%%% matrix format: 
%%% rows - stimulus events (e.g. location/identity); columns - responses 
datamat = round(100*rand(5, 5));

%%% to include all data in the fit
filtmat = ones(5, 5);

% %%% a "trick" to fit a forced-choice task (four alternative task) 
% %%% with catch trials, but *no* NoGo responses
% datamat(:,5) = 0.1;
% filtmat(:,5) = 0;

% %%% a "trick" to fit a forced-choice task (four alternative task) 
% %%% with *no* catch trials and *no* NoGo responses
% datamat(5,:) = 0.1;
% filtmat(5,:) = 0;

%%% to fit with only a subset of contingencies (e.g. last row/column) 
% filtmat = [0 0 0 0 1; ...
%            0 0 0 0 1; ...
%            0 0 0 0 1; ...
%            0 0 0 0 1; ...
%            1 1 1 1 1];
%
% % note, when including a subset of contingencies: 
% % 1) at least 8 independent contingencies must be included (8 params)
% % 2) at least one contingency in each row and column must be included
% % 3) at most four contingencies per row are independent (e.g., including 
% %    all five in a row amounts to only four independent observations)

%% perform ML estimation
disp('Optimization starts...')
[theta_est, theta_err, data_fit] = MLE_4ADC (datamat, filtmat);

%% post-processing

% sensitivity and bias parameters
sens = theta_est(1:4);
crit = theta_est(5:8);

% compute log-likelihood 
nobs     = sum(datamat(:));
prob_fit = data_fit/nobs;
LLFm     = datamat.*log(prob_fit);
LLF      = sum(LLFm(:));

%%% display
disp(' ');
disp('---------------------')
disp('Observed data')
disp(datamat)

disp('Fitted data')
disp(floor(data_fit))

disp(' ')
disp('ML estimated parameters')
fprintf('d1:%2.2f +/- %2.3f, d2:%2.2f +/- %2.3f, d3:%2.2f +/- %2.3f, d4:%2.2f +/- %2.3f \n', theta_est(1), theta_err(1), theta_est(2), theta_err(2), theta_est(3), theta_err(3), theta_est(4), theta_err(4)); 
fprintf('c1:%2.2f +/- %2.3f, c2:%2.2f +/- %2.3f, c3:%2.2f +/- %2.3f, c4:%2.2f +/- %2.3f \n'  , theta_est(5), theta_err(5), theta_est(6), theta_err(6), theta_est(7), theta_err(7), theta_est(8), theta_err(8)); 
disp('---------------------')
