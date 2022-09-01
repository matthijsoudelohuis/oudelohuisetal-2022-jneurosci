function [symbolstring] = MOL_BFtoSymbol(bf)
% Bayes Factor	Interpretation
% > 100	Extreme evidence for alternative hypothesis
% 30 – 100	Very strong evidence for alternative hypothesis
% 10 – 30	Strong evidence for alternative hypothesis
% 3 – 10	Moderate evidence for alternative hypothesis
% 1 – 3	Anecdotal evidence for alternative hypothesis
% 1	No evidence
% 1/3 – 1	Anecdotal evidence for null hypothesis
% 1/3 – 1/10	Moderate evidence for null hypothesis
% 1/10 – 1/30	Strong evidence for null hypothesis
% 1/30 – 1/100	Very strong evidence for null hypothesis
% < 1/100	Extreme evidence for null hypothesis
symbolstring = '';

if bf>=100 
    symbolstring = '****';
elseif bf>=30 && bf<100
    symbolstring = '***';
elseif bf>=10 && bf<30
    symbolstring = '**';
elseif bf>=3 && bf<10
    symbolstring = '*';
elseif bf>=1 && bf<3
    symbolstring = '';
elseif bf>=1/3 && bf<1
    symbolstring = '';
elseif bf>=1/10 && bf<1/3
    symbolstring = '#';
elseif bf>=1/30 && bf<1/10
    symbolstring = '##';
elseif bf>=1/100 && bf<1/30
    symbolstring = '###';
elseif bf<1/100
    symbolstring = '####';
end

end
