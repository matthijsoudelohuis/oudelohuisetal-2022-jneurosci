function [xvalsau,xvalsvis,ctable_fit_mat] = MOL_Gen2ADC_PsyCurve(theta,Par)

% Get curve from parameters:
% sivals = linspace(probepos,1,1000); %linear spacing
xvalsau = 10.^(linspace(log10(Par.auprobepos),log10(max(Par.auticks)),1000)); %logarithmic spacing
xvalsvis = 10.^(linspace(log10(Par.visprobepos),log10(max(Par.visticks)),1000)); %logarithmic spacing

dmax    = theta(1:2);
n       = theta(3:4);
s50     = theta(5:6);
ce      = theta(7:8);

% create a matrix of psychophysical functions for each location
dmat = nan(2, length(xvalsau));

% Compute d-prime for every value of visual and auditory:
dmat(1,:) = dmax(1) * (xvalsau.^n(1))./(xvalsau.^n(1) + s50(1)^n(1)); %Auditory
dmat(2,:) = dmax(2) * (xvalsvis.^n(2))./(xvalsvis.^n(2) + s50(2)^n(2)); %Visual

% for i = 1:2
%     dmat(i,:) = dmax(i) * (xvals.^n(i))./(xvals.^n(i) + s50(i)^n(i));
% end
ctable_fit_mat = NaN(3,3,length(xvalsau)); %init matrix
for i = 1:length(xvalsau) %change this to making directly from fit params
    ctable_fit_mat(:,:,i) = madc_pcont_table(dmat(:,i), ce, 'u');
end

end