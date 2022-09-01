function [ctable,auconditions,visconditions] = MOL_Prep_Psy2ADC(varargin)
%% Get input arguments:
tempsessionData     = varargin{1};
temptrialData       = varargin{2};

%% Get the trialtypes per session:
trialtypes = unique(temptrialData.trialType);

%% Define probe, auditory and visual trials and associated settings:
seshasprobe = 0; seshasau = 0; seshasvis = 0;
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
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasvis           = 1;
                visuallog          	= 1;
                vistrialtype        = 'V';
                vismultiply100      = 1;
                visxaxislabel       = 'Contrast (%)';
                vissessionfield     = 'vecVisStimIntensities';
                vistrialfield       = 'visualInt';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasvis           = 1;
                vismultiply100      = 0;
                visuallog           = 1;
                vistrialtype        = 'V';
                vissessionfield     = 'vecOriChange';
                vistrialfield       = 'visualOriChange';
                visxaxislabel       = 'Change in orientation (degrees)';
            end
        case 'X'
            seshasvis           = 1;
            vismultiply100      = 0;
            visuallog           = 1;
            vistrialtype        = 'X';
            vissessionfield     = 'vecOriChange';
            vistrialfield       = 'visualOriChange';
            visxaxislabel       = 'Change in orientation (degrees)';
        case 'A'
            if strcmp(tempsessionData.strTask,'mol_det') || strcmp(tempsessionData.strTask,'mol_mod_discr')
                seshasau            = 1;
                audiolog            = 0;
                aumultiply100       = 1;
                autrialtype         = 'A';
                ausessionfield      = 'vecAuStimIntensities';
                autrialfield        = 'audioInt';
                auxaxislabel        = 'Sound level (dB)';
            elseif strcmp(tempsessionData.strTask,'mol_ch_det')
                seshasau            = 1;
                aumultiply100       = 0;
                audiolog            = 1;
                autrialtype         = 'A';
                if strcmp(tempsessionData.auChangeUnit,'Hz')
                    ausessionfield      = 'vecFreqChange';
                    autrialfield        = 'audioFreqChange';
                    auxaxislabel        = 'Change in frequency (Hz)';
                else strcmp(tempsessionData.auChangeUnit,'Oct')
                    ausessionfield      = 'vecOctChange';
                    autrialfield        = 'audioOctChange';
                    auxaxislabel        = 'Change in octave (Oct)';
                end
            end
        case 'Y'
            seshasau            = 1;
            aumultiply100       = 0;
            audiolog            = 1;
            autrialtype         = 'Y';
            if strcmp(tempsessionData.auChangeUnit,'Hz')
                ausessionfield      = 'vecFreqChange';
                autrialfield        = 'audioFreqChange';
                auxaxislabel        = 'Change in frequency (Hz)';
            else strcmp(tempsessionData.auChangeUnit,'Oct')
                ausessionfield      = 'vecOctChange';
                autrialfield        = 'audioOctChange';
                auxaxislabel        = 'Change in octave (Oct)';
            end
    end
end

if seshasvis
    if iscell(tempsessionData.(vissessionfield))
        visconditions   = tempsessionData.(vissessionfield){1};
    else
        visconditions   = tempsessionData.(vissessionfield);
    end
end

if seshasau
    if iscell(tempsessionData.(ausessionfield))
        auconditions    = tempsessionData.(ausessionfield){1};
    else
        auconditions    = tempsessionData.(ausessionfield);
    end
end

%% Initialize output table:
if numel(visconditions)~=numel(auconditions)
    error('unequal number of conditions in modalities')
end
ctable = NaN(3,3,numel(visconditions));

%% Selection of responses per trial: (3 response options: visual/auditory/noresponse)
if tempsessionData.VisualLeftCorrectSide
    responseasauditory  = strcmp(temptrialData.responseSide,'R');
    responseasvisual    = strcmp(temptrialData.responseSide,'L');
    noresponse          = temptrialData.noResponse==1;
else
    responseasauditory  = strcmp(temptrialData.responseSide,'L');
    responseasvisual    = strcmp(temptrialData.responseSide,'R');
    noresponse          = temptrialData.noResponse==1;
end

%% Selection of probe trials:
selectionprobe      = strcmp(temptrialData.trialType,probetrialtype);
fprintf('%d Probe Trials\n',sum(selectionprobe))

%% Selection of auditory trials:
selectionau         = NaN(length(temptrialData.trialType),length(auconditions));
for i = 1:length(auconditions)
    selectionau(:,i)    = strcmp(temptrialData.trialType,autrialtype) & ismember(abs(temptrialData.(autrialfield)),auconditions(i));
    fprintf('%d Audio Trials (%d)\n',sum(selectionau(:,i)),auconditions(i))
end

%% %% Selection of visual trials:
selectionvis        = NaN(length(temptrialData.trialType),length(visconditions));
for i = 1:length(visconditions)
    selectionvis(:,i)    = strcmp(temptrialData.trialType,vistrialtype) & ismember(abs(temptrialData.(vistrialfield)),visconditions(i));
    fprintf('%d Visual Trials (%d)\n',sum(selectionvis(:,i)),visconditions(i))
end

%% Compute contingency table (intersection of stimulus/response per psy level):
for i = 1:length(auconditions)
    ctable(1,1,i)    = sum(selectionau(:,i) & responseasauditory) / sum(selectionau(:,i));
    ctable(1,2,i)    = sum(selectionau(:,i) & responseasvisual) / sum(selectionau(:,i));
    ctable(1,3,i)    = sum(selectionau(:,i) & noresponse) / sum(selectionau(:,i));
end

for i = 1:length(visconditions)
    ctable(2,1,i)    = sum(selectionvis(:,i) & responseasauditory) / sum(selectionvis(:,i));
    ctable(2,2,i)    = sum(selectionvis(:,i) & responseasvisual) / sum(selectionvis(:,i));
    ctable(2,3,i)    = sum(selectionvis(:,i) & noresponse) / sum(selectionvis(:,i));
end

ctable(3,1,:) = sum(selectionprobe & responseasauditory) / sum(selectionprobe);
ctable(3,2,:) = sum(selectionprobe & responseasvisual) / sum(selectionprobe);
ctable(3,3,:) = sum(selectionprobe & noresponse) / sum(selectionprobe);

end