function MOL_plotOptoBehaviorPsy_Rates_PPC(params,TotalResp)

figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.3 .5 .55]);    hold all;

% Settings:
params.auprobepos       = 1/64;
params.auticks          = [1/16 1/8 1/4 1/2];
% params.auticklabels     = {'Catch' '1/256' '1/64' '1/8' '1/2'};
params.auticklabels   = {'Catch' 'Sub' 'Thr' 'Sup' 'Max'};
params.auxaxislabel     = 'Change in octave';
params.auystatslabel    = 'Auditory threshold (partial octave)';

params.visprobepos      = 0.5;
params.visticks         = [2 8 30 90];
params.vistickslabels   = {'Catch' 'Sub' 'Thr' 'Sup' 'Max'};
params.visxaxislabel    = 'Change in orientation (Degrees)';
params.visystatslabel   = 'Visual threshold (Degrees)';

params.yticks           = [0 0.25 0.5 0.75 1];

% Construct x from position of the probe trial and conditions
xdata_vis               = [params.visprobepos params.visticks];

% Construct x from position of the probe trial and conditions
% xdata_au                = [params.auprobepos 10 50 260 4000];
xdata_au                = [params.auprobepos params.auticks];

% params.auticks         = [1/32 1/2];
% params.auticks         = [200 4000];
% params.auticklabels   = {'Probe' 'Thr' 'Max'};

% params.visticks         = [12 90];
% params.vistickslabels   = {'Probe' 'Thr' 'Max'};

%Audio:
subplot(1,2,1); hold all;
for iOpto = 1:3
    %Get the hit rates for audio conditions in this opto condition:
    %=full first dimension, vis=1 (no change), iOpto, resp=1,all sessions
    datatoplot = squeeze(TotalResp(:,1,iOpto,1,:));
%     plot(xdata_au,nanmean(datatoplot,2),'.','Color',params.colors_audio_opto{iOpto},'LineWidth',0.1,'MarkerSize',50);
    params.linehandles(iOpto) = errorbar(xdata_au,nanmean(datatoplot,2),nanstd(datatoplot,[],2)/sqrt(size(TotalResp,5)),params.lines_audio_opto{iOpto},'Color',params.colors_audio_opto{iOpto},...
        'MarkerFaceColor',params.colors_audio_opto{iOpto},'LineWidth',5,'MarkerSize',15);
end
legend(params.linehandles,{'Au - Control' 'Au - V1 Early' 'Au - V1 Late'},'FontSize',15,'Location','NorthEast');
legend boxoff;

for iCond=1:3
    datatotest          = squeeze(TotalResp(iCond,1,[1 2],1,:));
    tempbf              = bf.ttest(datatotest(1,:),datatotest(2,:));
    bfsymb              = MOL_BFtoSymbol(tempbf);
    text(xdata_au(iCond),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd               = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,2),tempd,tempbf)
end

for iOpto = 1:3
    %Get the hit rates for audio conditions in this opto condition:
    %=full first dimension, vis=1 (no change), iOpto, resp=1,all sessions
    datatoplot = squeeze(TotalResp(:,1,iOpto,2,:));
%     plot(xdata_au,nanmean(datatoplot,2),'.','Color',params.colors_audio_opto{iOpto},'LineWidth',0.1,'MarkerSize',50);
    params.linehandles(iOpto) = errorbar(xdata_au,nanmean(datatoplot,2),nanstd(datatoplot,[],2)/sqrt(size(TotalResp,5)),params.lines_visual_opto{iOpto},'Color',params.colors_visual_opto{iOpto},...
        'MarkerFaceColor',params.colors_visual_opto{iOpto},'LineWidth',5,'MarkerSize',15);
end
legend(params.linehandles,{'Ctrl' 'Opto'},'FontSize',15,'Location','NorthEast');
legend boxoff;

for iCond=1:3
    datatotest          = squeeze(TotalResp(iCond,1,[1 2],2,:));
    tempbf              = bf.ttest(datatotest(1,:),datatotest(2,:));
    bfsymb              = MOL_BFtoSymbol(tempbf);
    text(xdata_au(iCond),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd               = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,2),tempd,tempbf)
end


%%%%%%%%
%Visual:
subplot(1,2,2); hold all;
for iOpto = 1:3
    %Get the hit rates for visual conditions in this opto condition:
    %=au=1 (no change), full visual dimension, iOpto, resp=2,all sessions
    datatoplot = squeeze(TotalResp(1,:,iOpto,2,:));
    if size(datatoplot,1)==1; datatoplot = datatoplot'; end
    params.linehandles(iOpto) = errorbar(xdata_vis,nanmean(datatoplot,2),nanstd(datatoplot,[],2)/sqrt(size(TotalResp,5)),params.lines_visual_opto{iOpto},'Color',params.colors_visual_opto{iOpto},...
        'MarkerFaceColor',params.colors_visual_opto{iOpto},'LineWidth',5,'MarkerSize',15);
end
legend(params.linehandles,{'Ctrl' 'Opto'},'FontSize',15,'Location','NorthEast');
legend boxoff;

for iCond=1:3
    datatotest          = squeeze(TotalResp(1,iCond,[1 2],2,:));
    tempbf              = bf.ttest(datatotest(1,:),datatotest(2,:));
    bfsymb              = MOL_BFtoSymbol(tempbf);
    text(xdata_vis(iCond),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd               = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,2),tempd,tempbf)
end

for iOpto = 1:3
    %Get the hit rates for visual conditions in this opto condition:
    %=au=1 (no change), full visual dimension, iOpto, resp=2,all sessions
    datatoplot = squeeze(TotalResp(1,:,iOpto,1,:));
    if size(datatoplot,1)==1; datatoplot = datatoplot'; end
    params.linehandles(iOpto) = errorbar(xdata_vis,nanmean(datatoplot,2),nanstd(datatoplot,[],2)/sqrt(size(TotalResp,5)),params.lines_audio_opto{iOpto},'Color',params.colors_audio_opto{iOpto},...
        'MarkerFaceColor',params.colors_audio_opto{iOpto},'LineWidth',5,'MarkerSize',15);
end

for iCond=1:3
    datatotest          = squeeze(TotalResp(1,iCond,[1 2],1,:));
    tempbf              = bf.ttest(datatotest(1,:),datatotest(2,:));
    bfsymb              = MOL_BFtoSymbol(tempbf);
    text(xdata_vis(iCond),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd               = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,2),tempd,tempbf)
end

MOL_Psy2Sided_FigMakeup(params)


end