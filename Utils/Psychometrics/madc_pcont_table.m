function py = madc_pcont_table(d, c, mt)
%% generate contingency table response probabilities given d and c vector
%
%  Sridhar Devarajan, Sep 2016
%  sridhar@cns.iisc.ernet.in
%
% d and c should be mx1 each
% mt is model type: 'u' for unequal sensitivities and 'e' for equal sensitivities

if length(d) ~= length(c)
    error('d and c must be same length vectors');
end

m  = length(d);
py = nan(m+1, m+1);

global model_type
model_type = mt;

for i = 1:m, % m potential target locations
    
    %%% find response probabilities for each target contingency
    s = zeros(1,m);         % m potential target locations
    s(i) = 1; 
    
    % hits (i==r) and misidentifications (i~=r)
    for r = 1:m,
        Fe = @(x) pyr_integrand(x, c, d, r, s);
        py(i, r) = quadgk(Fe, c(r), inf);
    end
     
    % misses column
    py(i, m+1) = 1 - nansum(py(i, 1:m));
        
end

%  false alarm rate row
for r = 1:m
    s = zeros(1,m); 
    
    Fe = @(x) pyr_integrand(x, c, d, r, s);
    py(m+1, r) = quadgk(Fe, c(r), inf);
    
end

% correct rejection cell
py(m+1, m+1) = 1 - nansum(py(m+1, 1:m)); 


function pyr = pyr_integrand(x, c, d, r, s)
%% compute the madc integrand -- for any order of model
F = @normcdf; f = @normpdf;

pyr = ones(size(x)); m = length(s); 

global model_type

for i = 1:m,
    if i ~= r,
        if model_type == 'u'
            mv = F((d(r)/d(i))*(x - c(r)) + c(i), d(i)*s(i));   % unequal sensitivities
        elseif model_type == 'e'
            mv = F(x - c(r) + c(i), d(i)*s(i));                 % equal sensitivities
        end
    else
        mv = f(x, d(r)*s(r));
    end
    pyr = pyr .* mv;
end