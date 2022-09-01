function [theta_est, theta_err, data_fit] = MLE_4ADC (datamat, filtmat)
%% Maximum likelihood estimation of 4-ADC model parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datamat and filtmat are 5x5 matrices
% datamat - response counts (rows: stimulus events; columns: responses)
% filtmat - a switch for including/excluding contingencies in the fit
%
% (c) Sridharan Devarajan, Aug 2014
% dsridhar@stanford.edu; sridharan.d@gmail.com
%
% Citation: 
%   Sridharan, D., Steinmetz, N.A., Moore, T. and Knudsen, E.I. (2014).
%       Distinguishing bias from sensitivity effects in multialternative 
%       detection tasks. J. Vision 14 (9):16, 1-32. doi: 10.1167/14.9.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define probability densities and cumulative distributions (unit normal, by default) 
F1 = @normcdf; f1 = @normpdf;
F2 = F1; f2 = f1;
F3 = F1; f3 = f1;
F4 = F1; f4 = f1;

e_mu = 0; e_sig = 1; % mean, std of error


%% ----- nothing below this line needs editing ----- 
global modM
modM = {F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig};

global dmat         % 5x5 form
dmat = datamat;

global fmat         % 5x5 form
fmat = filtmat;

num_params = 8;
theta_init = rand(1, num_params);

[mLe_est, mLe_val, exitflag, output, grad, hessian] = fminunc(@evaluate_LLF, theta_init);

theta_est = mLe_est;
LLF       = -mLe_val; % -ve, because func value minimized

theta_se   = sqrt(inv(hessian));
theta_err  = diag(theta_se, 0); 


nX  = size(dmat, 1);
% ntr = repmat(sum(dmat,2)', nX, [])';
ntr = repmat(sum(dmat,2)', nX, 1)';

cX = [1];
prob_fit = get_data_fit(theta_est, cX);

prob_fit = reshape(prob_fit', [], nX)';
data_fit = prob_fit .* ntr; 

function F = evaluate_LLF(theta)
%% evaluate the log likelihood function

cX = [1];

Xeye = eye(4);
X1 = Xeye(1,:);
X2 = Xeye(2,:);
X3 = Xeye(3,:);
X4 = Xeye(4,:);
X0 = zeros(1,4);

global modM

[p1] = generate_parametric_curves(cX, X1, theta, modM);
[p2] = generate_parametric_curves(cX, X2, theta, modM);
[p3] = generate_parametric_curves(cX, X3, theta, modM);
[p4] = generate_parametric_curves(cX, X4, theta, modM);
[p0] = generate_parametric_curves(cX, X0, theta, modM);

[p11, p21, p31, p41, p01] = deal(p1{:});
[p12, p22, p32, p42, p02] = deal(p2{:});
[p13, p23, p33, p43, p03] = deal(p3{:});
[p14, p24, p34, p44, p04] = deal(p4{:});
[p10, p20, p30, p40, p00] = deal(p0{:});

pij = [p11, p21, p31, p41, p01, p12, p22, p32, p42, p02, p13, p23, ...
    p33, p43, p03, p14, p24, p34, p44, p04, p10, p20, p30, p40, p00]; % here prs => p_s^r (prob resp to r when stim at s)

global dmat
global fmat

yij = dmat;

nX = 5;
pij = reshape(pij', [], nX)';
LLF = 0;

for i = 1:nX
    y0ij = yij(i, fmat(i,:) == 0);
    y1ij = yij(i, fmat(i,:) == 1);
    
    p0ij = pij(i, fmat(i,:) == 0);
    p1ij = pij(i, fmat(i,:) == 1);
    
    LLF = LLF + sum(y1ij.*log(p1ij));
    
    if ~isempty(p0ij)
        LLF = LLF + sum(y0ij)*log(sum(p0ij)); 
    end
end

F   = -LLF; % minimize negative likelihood => maximize likelihood
            % note: this is not the actual LLF, it misses some constant
            % (additive) terms, but will do for our purposes


function data_fit = get_data_fit(theta, cX)
%% evaluate probabilities with model parameters

Xeye = eye(4);
X1 = Xeye(1,:);
X2 = Xeye(2,:);
X3 = Xeye(3,:);
X4 = Xeye(4,:);
X0 = zeros(1,4);

global modM

[p1] = generate_parametric_curves(cX, X1, theta, modM);
[p2] = generate_parametric_curves(cX, X2, theta, modM);
[p3] = generate_parametric_curves(cX, X3, theta, modM);
[p4] = generate_parametric_curves(cX, X4, theta, modM);
[p0] = generate_parametric_curves(cX, X0, theta, modM);

[p11, p21, p31, p41, p01] = deal(p1{:});
[p12, p22, p32, p42, p02] = deal(p2{:});
[p13, p23, p33, p43, p03] = deal(p3{:});
[p14, p24, p34, p44, p04] = deal(p4{:});
[p10, p20, p30, p40, p00] = deal(p0{:});

data_fit = [p11, p21, p31, p41, p01, p12, p22, p32, p42, p02, p13, p23, p33, p43, p03, p14, p24, p34, p44, p04, p10, p20, p30, p40, p00]; 


%% internal functions (that evaluate the integral)
function [x1, x2, x3, x4, ...
    d1, d2, d3, d4, c1, c2, c3, c4, ...
        F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] = dealargs(X, theta, modM)
xC = num2cell(X);
thetaC = num2cell(theta);

[x1, x2, x3, x4] = deal(xC{:});
[d1, d2, d3, d4, c1, c2, c3, c4] = deal(thetaC{:});
[F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] = deal(modM{:});

function py1 = comput_integ_y1(X, theta, modM)
[x1, x2, x3, x4, ...
    d1, d2, d3, d4, c1, c2, c3, c4, ...
        F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] ...
            = dealargs(X, theta, modM);

Llim = c1 - d1*x1; Ulim = Inf; 
F_e1 = @(x)F2(x + d1*x1 - c1 - d2*x2 + c2, e_mu, e_sig).*F3(x + d1*x1 - c1 - d3*x3 + c3, e_mu, e_sig) ...
            .*F4(x + d1*x1 - c1 - d4*x4 + c4, e_mu, e_sig).*f1(x, e_mu, e_sig);
py1 = quadgk(F_e1, Llim, Ulim);


function py2 = comput_integ_y2(X, theta, modM)
[x1, x2, x3, x4, ...
    d1, d2, d3, d4, c1, c2, c3, c4, ...
        F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] ...
            = dealargs(X, theta, modM);

Llim = c2 - d2*x2; Ulim = Inf; 
F_e2 = @(x)F1(x + d2*x2 - c2 - d1*x1 + c1, e_mu, e_sig).*F3(x + d2*x2 - c2 - d3*x3 + c3, e_mu, e_sig) ...
            .*F4(x + d2*x2 - c2 - d4*x4 + c4, e_mu, e_sig).*f2(x, e_mu, e_sig);
py2 = quadgk(F_e2, Llim, Ulim);

function py3 = comput_integ_y3(X, theta, modM)
[x1, x2, x3, x4, ...
    d1, d2, d3, d4, c1, c2, c3, c4, ...
        F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] ...
            = dealargs(X, theta, modM);

Llim = c3 - d3*x3; Ulim = Inf; 
F_e3 = @(x)F1(x + d3*x3 - c3 - d1*x1 + c1, e_mu, e_sig).*F2(x + d3*x3 - c3 - d2*x2 + c2, e_mu, e_sig) ...
            .*F4(x + d3*x3 - c3 - d4*x4 + c4, e_mu, e_sig).*f3(x, e_mu, e_sig);
py3 = quadgk(F_e3, Llim, Ulim);

function py4 = comput_integ_y4(X, theta, modM)
[x1, x2, x3, x4, ...
    d1, d2, d3, d4, c1, c2, c3, c4, ...
        F1, F2, F3, F4, f1, f2, f3, f4, e_mu, e_sig] ...
            = dealargs(X, theta, modM);

Llim = c4 - d4*x4; Ulim = Inf; 
F_e4 = @(x)F1(x + d4*x4 - c4 - d1*x1 + c1, e_mu, e_sig).*F2(x + d4*x4 - c4 - d2*x2 + c2, e_mu, e_sig) ...
            .*F3(x + d4*x4 - c4 - d3*x3 + c3, e_mu, e_sig).*f4(x, e_mu, e_sig);
py4 = quadgk(F_e4, Llim, Ulim);

function [p_y_x] = generate_parametric_curves(cX, X, theta, modM)

p_y1_x = nan(1, length(cX));
p_y2_x = nan(1, length(cX));
p_y3_x = nan(1, length(cX));
p_y4_x = nan(1, length(cX));
p_y0_x = nan(1, length(cX));

for i_cX = 1:length(cX),
    p_y1_x(i_cX) = comput_integ_y1(cX(i_cX)*X, theta, modM);
    p_y2_x(i_cX) = comput_integ_y2(cX(i_cX)*X, theta, modM);
    p_y3_x(i_cX) = comput_integ_y3(cX(i_cX)*X, theta, modM);
    p_y4_x(i_cX) = comput_integ_y4(cX(i_cX)*X, theta, modM);
    p_y0_x(i_cX) = 1 - p_y1_x(i_cX) - p_y2_x(i_cX) - p_y3_x(i_cX) - p_y4_x(i_cX);
end

p_y_x = {p_y1_x, p_y2_x, p_y3_x, p_y4_x, p_y0_x};