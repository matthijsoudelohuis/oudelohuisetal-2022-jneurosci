function [visconditions,auconditions,FullRespMat,FullnTrialsMat] = MOL_Psy_GetTrialMat(varargin)
%Take trialdata and sessiondata as input and creates a stimulus x response outcome matrix:
% Output is a M (rows) by N (columns) by P matrix,
% where M is the number of auditory conditions (row 1 is no auditory change)
% where N is the number of visual conditions (column 1 is no visual change)
% and P is the number of response options (1 is auditory, 2 is visual, 3 is no response)

%% Get input arguments:
tempsessionData     = varargin{1};
temptrialData       = varargin{2};

%% Get the trialtypes per session:
trialtypes = unique(temptrialData.trialType);

%% Define probe, auditory and visual trials and associated settings:
seshasprobe = 0; seshasau = 0; seshasvis = 0; seshasconflict = 0;
for tt = 1:length(trialtypes)
    switch trialtypes{tt}
        case 'P'
            seshasprobe         = 1;
            probetrialtype      = 'P';
        case 'Q'
            seshasprobe         = 1;
            probetrialtype      = 'Q';
        case 'R'
            seshasprobe         = 1;
            probetrialtype      = 'R';
        case 'V'
            if all(strcmp(tempsessionData.strTask,'mol_det') | strcmp(tempsessionData.strTask,'mol_mod_discr'))
                seshasvis           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecVisStimIntensities';
                vistrialfield       = 'visualInt';
            elseif all(strcmp(tempsessionData.strTask,'mol_ch_det'))
                seshasvis           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecOriChange';
                vistrialfield       = 'visualOriChange';
            end
        case 'X'
            seshasvis           = 1;
            vistrialtype        = 'X';
            vissessionfield     = 'vecOriChange';
            vistrialfield       = 'visualOriChange';
        case 'A'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasau            = 1;
                autrialtype         = 'A';
                ausessionfield      = 'vecAuStimIntensities';
                autrialfield        = 'audioInt';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasau            = 1;
                autrialtype         = 'A';
                if strcmp(tempsessionData.auChangeUnit,'Hz') || isempty(tempsessionData.auChangeUnit{1})
                    ausessionfield      = 'vecFreqChange';
                    autrialfield        = 'audioFreqChange';
                else strcmp(tempsessionData.auChangeUnit,'Oct')
                    ausessionfield      = 'vecOctChange';
                    autrialfield        = 'audioOctChange';
                end
            end
        case 'Y'
            seshasau            = 1;
            autrialtype         = 'Y';
            if strcmp(tempsessionData.auChangeUnit(1),'Hz') || isempty(tempsessionData.auChangeUnit{1})
                ausessionfield      = 'vecFreqChange';
                autrialfield        = 'audioFreqChange';
            elseif strcmp(tempsessionData.auChangeUnit(1),'Oct')
                ausessionfield      = 'vecOctChange';
                autrialfield        = 'audioOctChange';
            end
        case 'C'
            seshasconflict      = 1;
            conflicttrialtype   = 'C';
    end
end

%% Visual:
if seshasvis
    if iscell(tempsessionData.(vissessionfield))
        visconditions   = tempsessionData.(vissessionfield){1};
    else
        visconditions   = tempsessionData.(vissessionfield);
    end
    visconditions = visconditions(visconditions~=0);  visconditions = visconditions(:)';
    
end

%% Auditory conditions
if seshasau
    if iscell(tempsessionData.(ausessionfield))
        auconditions    = tempsessionData.(ausessionfield){1};
    elseif size(tempsessionData.(ausessionfield),1)>1
        auconditions    = tempsessionData.(ausessionfield)(1,:);
    else
        auconditions    = tempsessionData.(ausessionfield);
    end
    auconditions = auconditions(auconditions~=0); auconditions = auconditions(:)';
    
end

%% Init matrices to store output:
FullRespMat     = NaN(length(auconditions)+1,length(visconditions)+1,3);
FullnTrialsMat  = NaN(length(auconditions)+1,length(visconditions)+1);

%% Selection of probe trials:
if seshasprobe
    selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
else selectionprobe     = strcmp(temptrialData.trialType,'Z'); %Nonsensical index
end

%% Get responses:
if seshasprobe
    FullnTrialsMat(1,1)         = sum(selectionprobe);
    for iResp = 1:3
        FullRespMat(1,1,iResp)          = sum(temptrialData.vecResponse==iResp & selectionprobe);
    end
end

if seshasau % Auditory only trials:
    for iAu = 1:length(auconditions)
        idx = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(iAu));
        FullnTrialsMat(iAu+1,1)    = sum(idx);
        for iResp = 1:3
            FullRespMat(iAu+1,1,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
        end
    end
end

if seshasvis % Visual only trials:
    for iVis = 1:length(visconditions)
        idx = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(iVis));
        FullnTrialsMat(1,iVis+1)    = sum(idx);
        for iResp = 1:3
            FullRespMat(1,iVis+1,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
        end
    end
end

if seshasconflict
    for iAu = 1:length(auconditions)
        for iVis = 1:length(visconditions)
            %Get index of trials as crosssection of both au and vis conditions:
            idx = ismember(abs(temptrialData.(autrialfield)),auconditions(iAu))...
                & ismember(abs(temptrialData.(vistrialfield)),visconditions(iVis))...
                & strcmp(temptrialData.trialType,conflicttrialtype);
            FullnTrialsMat(iAu+1,iVis+1)    = sum(idx);
            for iResp = 1:3
                FullRespMat(iAu+1,iVis+1,iResp)          = sum(temptrialData.vecResponse==iResp & idx);
            end
        end
    end
end

%% Checks and balances:
if FullnTrialsMat ~= sum(FullRespMat,3)
    error('trials dont match')
end

FullRespMat             = FullRespMat ./ repmat(FullnTrialsMat,1,1,3); %Make fractional (as percentage correct etc.)

end