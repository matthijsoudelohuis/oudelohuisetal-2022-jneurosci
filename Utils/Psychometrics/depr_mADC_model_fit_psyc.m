function [theta_est, theta_err, LLF, ctable_mod, ctable_fit, sivals, psyc_mate, ce] = mADC_model_fit_psyc(ctable, cvals, d_limit)

%% function that fits the m-ADC model to data in a contingency table
% rather than fitting d's to individual stimulus strengths, fits the entire
% psychophysical function in one shot
%
%  Sridhar Devarajan, Sep 2016
%  sridhar@cns.iisc.ernet.in
%
% Usage:
%   [theta_est, theta_err, LLF, ctable_fit] = mADC_model_fit_psyc(ctable, cvals)
%
% ctable is a 3D contingency table of the form (m+1) x (m+1) x k, 
% cvals is a 1D array of stimulus strength of size 1xk; *must* be normalized
% from 0-1, where 1 is the expected stimulus strength for asymptotic dmax
%
% where: 
%   m is the number of stimulus/response locations + 1 catch/NoGo
%   k is the number of stimulus strengths tested (e.g. contrast)
%   
%   data are expected to be arranged in increasing order of stimulus strength
%   the last row of each ctable should be identical and contain the false
%   alarm counts across all locations


% model type: 'u' for unequal sensitivities and 'e' for equal sensitivities
global model_type
model_type = 'u';   

if size(ctable, 1) ~= size(ctable, 2) 
    error('number of rows and columns of contingency table have to match');
end

if size(ctable, 3) ~= length(cvals)
    error('third dimension of ctable should be same length as cvals');
end
    
global M K
M = size(ctable, 1) - 1;
K = size(ctable, 3);

global yobs

yobs= [];

for i = 1:K
    yobs{i} = squeeze(ctable(:,:,i));
    yobs{i}(end,:) = (1/K) * yobs{i}(end,:);        % we divide the false alarm counts (last row of each table) by the number of stimulus strengths -- this prevents "overfitting" to the false alarm rates
end

global svals
svals = cvals;

%%% start the fitting

%% initialize
num_params = M*4;   % 3 d' -related parameters and one 'c' parameter for each location
                    % 3 d' parameters are dmax, n and s50, for each location
                    
dmax_init = 2.0; n_init = 1.0; s50_init = median(svals);
theta_init = [dmax_init*ones(1,M), n_init*ones(1,M), s50_init*ones(1,M) rand(1, M)];

%% constraints
if ~isempty(d_limit)
    dmax_lb = d_limit(1); dmax_ub = d_limit(2);
else
    dmax_lb = 1e-4;  dmax_ub = 5.0;
end

n_lb    = 1.0;  n_ub    = 3.0;
s50_lb  = 0.0;  s50_ub  = max(svals);

c_lb = -5.0;    c_ub = 5.0;

lb = [dmax_lb * ones(1,M), n_lb * ones(1,M), s50_lb * ones(1,M), c_lb * ones(1,M)];
ub = [dmax_ub * ones(1,M), n_ub * ones(1,M), s50_ub * ones(1,M), c_ub * ones(1,M)];

%%% additional  constraints here
%Aeq = [zeros(1,m*k) zeros(1,k)]; Beq = [0.0];

%% fitting
[mLe_est, mLe_val, exitflag, output, lambda, grad, hessian] = fmincon(@madc_LLF, theta_init, [], [], [], [], lb, ub);

theta_est = mLe_est;
LLF       = -mLe_val; % -ve, because func value minimized, so LLF is maximized

theta_se   = sqrt(inv(hessian));
theta_err = diag(theta_se, 0);

%% the best estimate and parameters
py = get_madc_probs(theta_est);

ctable_fit = cell(1,K);
for i=1:K
    ctable_fit{i} = py{i} .* repmat(sum(yobs{i}, 2), 1, M+1);
end

ctable_mod = yobs;

%% let's try some plotting

dmaxe   = theta_est(1:M); 
ne      = theta_est(M+1:2*M); 
s50e    = theta_est(2*M+1:3*M);
ce      = theta_est(3*M+1:4*M);

[psyc_mate, sivals] = get_interp_psycmat(dmaxe, ne, s50e);
% 
% figure(); set(gcf, 'Position', [100 100 350 300], 'Color', [1 1 1]); ylim([0,5]); hold on;
% colc = hsv(M);
% 
% for i=1:M
%     plot(sivals, psyc_mate(i,:), 'Color', colc(i,:), 'Linewidth', 2); 
%     plot([min(sivals), max(sivals)], [ce(i), ce(i)], 'Linestyle', '--', 'Color', colc(i,:), 'Linewidth', 1); 
%     box off;
% end
% 
% set(gca, 'Fontsize', 18); xlabel('\Delta Theta'); ylabel('d'' or c');



function F = madc_LLF (theta)
%% evaluate the LLF based on fit of real data to cont table produced by these values of d/c

global yobs K

py = get_madc_probs(theta);

LLF = 0;
for i = 1:K,    
    LLF = LLF + sum(sum(yobs{i} .* log(py{i})));
end

F   = -LLF;

function py = get_madc_probs(theta, int_flag)
% return a k element cell array with an m+1 x m+1 contingency table for
% each stimulus strength

global model_type M K

dmax = theta(1:M); n = theta(M+1:2*M); s50 = theta(2*M+1:3*M);

psyc_mat = get_psycmat(dmax, n, s50);

py = cell(1,K);
for i = 1:K,
    d = psyc_mat(:, i)';  % get the d' value from the psychophysical function matrix for that stimulus strength (i) for all locations 1:M
    c = theta(end-M+1:end);

    py{i} = madc_pcont_table(d, c, model_type);    
end

function dmat = get_psycmat(dmax, n, s50)
% create an matrix of psychophysical functions for each location
global M K svals

dmat = nan(M, K);
for i = 1:M
    dmat(i,:) = dmax(i) * (svals.^n(i))./(svals.^n(i) + s50(i)^n(i));
end

function [dmat, sivals] = get_interp_psycmat(dmax, n, s50)
% create an matrix of psychophysical functions for each location
global M K svals

sivals = linspace(min(svals), max(svals), 1000);
dmat = nan(M, length(sivals));
for i = 1:M
    dmat(i,:) = dmax(i) * (sivals.^n(i))./(sivals.^n(i) + s50(i)^n(i));
end