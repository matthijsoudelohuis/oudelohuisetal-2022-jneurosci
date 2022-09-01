function [params] = MOL_getColors_CHDET(params)

%% 
params.colors_modalities    = {[1 0 0] [0 0 1]};        %Auditory and Visual resp.

params.colors_trialtypes    = {[0.8 0 0] [0 0 0.8] [0.3 0.3 0.3] [0.8 0 0.8]}; %A, V, P, C

params.colors_conflict_choices = {[224 11 120] [95 15 148]}; %AC-A, C-V
params.colors_conflict_choices = cellfun(@(x) x/256,params.colors_conflict_choices,'UniformOutput',false);

%% Different trial types with opto:
params.colors_visual_opto   = {[0 0 0.9] [0.6 0.6 1] [0.4 0.4 1]};
params.colors_audio_opto    = {[0.9 0 0] [1 0.6 0.6] [1 0.4 0.4]};

%%
params.colors_ztrials       = {[7 128 148] [33 7 138] [148 134 15] [148 74 15]}; 
params.colors_ztrials       = cellfun(@(x) x/256,params.colors_ztrials,'UniformOutput',false);

%% Different experiments/cohorts/tasks:
params.colors_experiments   = {[19 160 67] [253 89 27] [116 25 114]}; params.colors_experiments = cellfun(@(x) x/256,params.colors_experiments,'UniformOutput',false);

%% %Dprime and criterion
params.colors_params        = {[0 0 0] [0.4 0.2 0.3]};  %Dprime and criterion

%% State of the animal
params.colors_state         = {[0 168 0] [168 0 84]};  %Dprime and criterion
params.colors_state         = cellfun(@(x) x/256,params.colors_state,'UniformOutput',false);
params.labels_state         = {'Passive' 'Active'};

%% AUC

params.labels_AUC            = {'Ori-small' 'Ori-large' 'Freq-small' 'Freq-large'...
    'Vis-small' 'Vis-large' 'Aud-small' 'Aud-large'...
    'Vis-small-outcome' 'Vis-large-outcome' 'Aud-small-outcome' 'Aud-large-outcome'...
    'Resp-sm12'  'Resp-la12' 'Resp-sm34' 'Resp-la34'...
    'Resp-sm12'  'Resp-la12' 'Resp-sm34' 'Resp-la34'...
    'Ori-small-hit' 'Ori-small-miss' 'Ori-large-hit' 'Ori-large-miss' 'Ori-hit' 'Ori-miss'...
    'Vis-small' 'Vis-large' 'Aud-small' 'Aud-large'...
    'Vis-all' 'Aud-all' 'Hit-probe'};
params.colors_AUC            = {[1 0 1]      [0.6 0 0.6]     [1 0 1]     [0.6 0 0.6]...
    [0.3 0.3 1]     [0 0 1]       [1 0.3 0.3]     [1 0 0]...
    [0.5 0.5 0.5]   [0 0 0]         [0.5 0.5 0.5] [0 0 0]...
    [0.5 0.5 0.5]   [0 0 0]         [0.5 0.5 0.5] [0 0 0] ...
    [0.5 0.5 0.5]   [0 0 0]         [0.5 0.5 0.5] [0 0 0] ...
    [1 0 1]         [0.6 0 0.6]     [1 0 1]     [0.6 0 0.6]  [1 0 1] [0.6 0 0.6] ...
    [0.3 0.3 1]     [0 0 1]         [1 0.3 0.3]     [1 0 0]...
    [0 0 1] [1 0 0] [0.2 0.2 0.2]};


end