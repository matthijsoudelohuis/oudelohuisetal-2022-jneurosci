function MOL_Psy2Sided_FigMakeup(Par)

%% Auditory subplot:
subplot(1,2,1); hold all;
xlim([Par.auprobepos max(Par.auticks)*1.02])
ylim([0 1])
xlabel(Par.auxaxislabel,'FontSize', 20)
ylabel('Hit Rate (%)','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','left','linewidth',2)
set(gca,'XScale','log');
set(gca,'Xdir','reverse','XTick',[Par.auprobepos Par.auticks],'XTickLabels',Par.auticklabels);
set(gca,'YTick',Par.yticks,'YTickLabels',Par.yticks*100);
set(gca,'XMinorTick','off')
box on
if isfield(Par,'aulinehandles')
    legend(Par.aulinehandles,{'Auditory response' 'Visual response'},'FontSize',15,'Location','NorthEast');
    legend(gca,'boxoff');
end

%% Visual Subplot:
subplot(1,2,2); hold all;
xlim([Par.visprobepos max(Par.visticks)*1.02])
ylim([0 1])
xlabel(Par.visxaxislabel,'FontSize', 20)
ylabel('Hit Rate (%)','FontSize', 20)
set(gca,'FontSize',15,'YAxisLocation','right')
set(gca,'linewidth',2)
set(gca,'XScale','log');
set(gca,'YTick',Par.yticks,'YTickLabels',Par.yticks*100);
set(gca,'Xdir','normal','XTick',[Par.visprobepos Par.visticks],'XTickLabels',Par.vistickslabels);
set(gca,'XMinorTick','off')
box on
if isfield(Par,'vislinehandles')
    legend(Par.vislinehandles,{'Visual response' 'Auditory response'},'FontSize',15,'Location','NorthWest');
    legend(gca,'boxoff');
end

end