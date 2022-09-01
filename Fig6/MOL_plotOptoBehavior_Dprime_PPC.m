function MOL_plotOptoBehavior_Dprime_PPC(params,dVis,dAud)
%% Dprime figure:
% Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.1 0.3 .3 .57]); hold all;
y_mean      = squeeze(nanmean(dVis(:,1,:))); y_mean = y_mean(~isnan(y_mean));
y_std       = squeeze(nanstd(dVis(:,1,:))) / sqrt(sum(~isnan(dVis(:,1,1)))); y_std = y_std(~isnan(y_std));
xpos        = 1:length(y_mean);
audoffset   = length(y_mean)+1;
errorbar(xpos, y_mean,y_std,':ob','MarkerSize',20,'MarkerEdgeColor','blue','LineWidth',3);

y_mean      = squeeze(nanmean(dVis(:,2,:)));  y_mean = y_mean(~isnan(y_mean));
y_std       = squeeze(nanstd(dVis(:,2,:))) / sqrt(sum(~isnan(dVis(:,2,1)))); y_std = y_std(~isnan(y_std));
errorbar(xpos, y_mean,y_std,'-ob','MarkerSize',20,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',3);

y_mean      = squeeze(nanmean(dAud(:,1,:))); y_mean = y_mean(~isnan(y_mean));
y_std       = squeeze(nanstd(dAud(:,1,:))) / sqrt(sum(~isnan(dAud(:,1,1))));  y_std = y_std(~isnan(y_std));
errorbar(xpos + audoffset, y_mean,y_std,':or','MarkerSize',20,'MarkerEdgeColor','red','LineWidth',3)

y_mean      = squeeze(nanmean(dAud(:,2,:))); y_mean = y_mean(~isnan(y_mean));
y_std       = squeeze(nanstd(dAud(:,2,:))) / sqrt(sum(~isnan(dAud(:,2,1)))); y_std = y_std(~isnan(y_std));
errorbar(xpos + audoffset, y_mean,y_std,'-or','MarkerSize',20,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',3)

barlocations = 1:2;
for iTrial = 1:2
    datatotest          = squeeze(dVis(:,iTrial,:));
    %Statistical testing:   wilcoxon signed rank test:
%     p = signrank(datatotest(:,1),datatotest(:,2));
%     sigstar(barlocations,p) %use sigstar function to identify sign conditions
    
    tempbf      = bf.ttest(datatotest(:,1),datatotest(:,2));
    bfsymb      = MOL_BFtoSymbol(tempbf);
    text(mean(barlocations),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd       = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,1),tempd,tempbf)
%      (i.e. the number of shocks in the coerced condition, t(36)=?3.177, p = 0.003, Cohen's d=?0.522, BF10=11.70). 
     
end

barlocations = 4:5;
for iTrial = 1:2
    datatotest          = squeeze(dAud(:,iTrial,:));
    %Statistical testing:   wilcoxon signed rank test:
%     p = signrank(datatotest(:,1),datatotest(:,2));
%     sigstar(barlocations,p) %use sigstar function to identify sign conditions
    
    tempbf      = bf.ttest(datatotest(:,1),datatotest(:,2));
    bfsymb      = MOL_BFtoSymbol(tempbf);
    text(mean(barlocations),nanmean(datatotest(:)),bfsymb,'FontSize',15)
    tempd       = computeCohen_d(datatotest(:,1),datatotest(:,2),'paired');
    fprintf('Paired Bayesian ttest: %d sessions, Cohen''s d = %1.3f, BF10=%3.2f\n',size(datatotest,1),tempd,tempbf)
end

%Make up:
ylim([-0.2 3])
ylabel('Dprime')
XTickLabels = repmat({'Control' 'Full Inh'},1,2);
set(gca,'XTick',[xpos xpos+audoffset],'XTickLabels',XTickLabels,'XTickLabelRotation',60);
grid on;
legend({'Small visual change' 'Large visual change' 'Small auditory change' 'Large auditory change'},'Location','north');
legend boxoff

end