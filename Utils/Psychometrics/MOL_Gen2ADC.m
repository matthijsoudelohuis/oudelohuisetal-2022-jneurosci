function outcome = MOL_Gen2ADC(d1,d2,c1,c2)

%% Generative Model:
% d1 = dprime signal 1
% d2 = dprime signal 2
% c1 = criterion signal 1
% c2 = criterion signal 2

X1stim = [1 0 0];
X2stim = [0 1 0];

outcome = NaN(3,3);

for signal = 1:3
    %Set signal for the three possible conditions:
    X1 = X1stim(signal);
    X2 = X2stim(signal);
    
    f1 = @(x) (normcdf(x+d1*X1-d2*X2-c1+c2)) .* (exp(-0.5 * ((x - 0)./1).^2) ./ (sqrt(2*pi) .* 1)); %Noise distribution
    outcome(signal,1)       = integral(f1,c1-d1*X1,Inf);
    
    f2 = @(x) (normcdf(x+d2*X2-d1*X1-c2+c1)) .* (exp(-0.5 * ((x - 0)./1).^2) ./ (sqrt(2*pi) .* 1)); %Noise distribution
    outcome(signal,2)       = integral(f2,c2-d2*X2,Inf);
   
    outcome(signal,3)       = normcdf(c1-d1*X1,0,1) * normcdf(c2-d2*X2,0,1);

end

% Show the discrete outcome in response probabilities for the 3X3 matrix:
disp(' ');
disp('           |              Response                 |')
disp('  Signal   |  "1"         |  "2"       |  "Catch"')
disp('  --------------------------------------------------|')
fprintf('  "1"      |   %3.1f      |   %3.1f       |   %3.1f    | \n',100*outcome(1,1),100*outcome(1,2),100*outcome(1,3));
disp('  ---------+--------------+-------------------------|');
fprintf('  "2"      |   %3.1f      |   %3.1f       |   %3.1f    |\n',100*outcome(2,1),100*outcome(2,2),100*outcome(2,3));
disp('  --------------------------------------------------|');
fprintf('  Absent   |   %3.1f      |   %3.1f       |   %3.1f    |\n',100*outcome(3,1),100*outcome(3,2),100*outcome(3,3));
disp('  --------------------------------------------------');


end