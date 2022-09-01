function [theta_est, theta_err, LLF, ctable_mod, ctable_fit] = mADC_model_fit(ctable, d_limit)

%% function that fits the m-ADC model to data in a contingency table
%
%  Sridhar Devarajan, Sep 2016
%  sridhar@cns.iisc.ernet.in
%
% Usage:
%   [theta_est, theta_err, LLF, ctable_fit] = mADC_model_fit(ctable)
%
% ctable is a 3D contingency table of the form (m+1) x (m+1) x k,
% where:
%   m is the number of stimulus/response locations + 1 catch/NoGo
%   k is the number of stimulus strengths tested (e.g. contrast)
%   data are expected to be arranged in increasing order of stimulus strength
%   the last row of each ctable should be identical and contain the false
%   alarm counts across all locations
%

% model type: 'u' for unequal sensitivities and 'e' for equal sensitivities
global model_type
model_type = 'u';

if size(ctable, 1) ~= size(ctable, 2)
    error('number of rows and columns of contingency table have to match');
end
global M K
M = size(ctable, 1) - 1;
K = size(ctable, 3);

global yobs

yobs = [];
for i = 1:K
    yobs{i} = squeeze(ctable(:,:,i));
    yobs{i}(end,:) = (1/K) * yobs{i}(end,:);         % we divide the false alarm counts (last row of each table) by the number of stimulus strengths -- this prevents "overfitting" to the false alarm rates
end


%%% start the fitting

%% initialize
num_params = M*(K+1);   % 'k' d parameters and one 'c' parameter for each location
theta_init = rand (1, num_params);

%% constraints

if nargin > 1    
    if ~isempty(d_limit)
        d_lb = d_limit(1);    d_ub = d_limit(2);
    else
        d_lb = 1e-4;    d_ub = 5.0;
    end
else
    d_lb = 1e-4;    d_ub = 5.0;
end

c_lb = -5.0;    c_ub = 5.0;

lb = [d_lb * ones(1,M*K), c_lb * ones(1,M)];
ub = [d_ub * ones(1,M*K), c_ub * ones(1,M)];

%%% additional  constraints here
%Aeq = [zeros(1,m*k) zeros(1,k)]; Beq = [0.0];

%% fitting
%options = optimoptions('fmincon','MaxFunctionEvaluations',1000,'MaxIterations',3000,'Display','iter');
[mLe_est, mLe_val, exitflag, output, lambda, grad, hessian] = fmincon(@madc_LLF, theta_init, [], [], [], [], lb, ub,[],[]);

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

function F = madc_LLF (theta)
%% evaluate the LLF based on fit of real data to cont table produced by these values of d/c

global yobs K

py = get_madc_probs(theta);

LLF = 0;
for i = 1:K,
    LLF = LLF + sum(sum(yobs{i} .* log(py{i})));
end

F   = -LLF;

function py = get_madc_probs(theta)
% return a k element cell array with an m+1 x m+1 contingency table for
% each stimulus strength

global model_type M K

py = cell(1,K);

for i = 1:K,
    d = theta((i-1)*M+1:i*M);
    c = theta(end-M+1:end);
    
    py{i} = madc_pcont_table(d, c, model_type);
end


%% Old code
% % % if ~exist('c_init', 'var') && ~exist('d_init', 'var')
% % %     theta_init = rand(1, 2*m);
% % % elseif ~exist('c_init', 'var') && exist('d_init', 'var')
% % %     theta_init = [d_init, rand(1, m)];
% % % elseif exist('c_init', 'var') && ~exist('d_init', 'var')
% % %     theta_init = [rand(1, m), c_init];
% % % end
