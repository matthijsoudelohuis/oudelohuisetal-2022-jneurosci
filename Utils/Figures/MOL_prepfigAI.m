function MOL_prepfigAI()

allAxesInFigure = findall(gcf,'type','axes');

for i = 1:length(allAxesInFigure)
    pause(0.1)
    ax = allAxesInFigure(i);
    axes(ax)
    xlims = xlim();
    ylims = ylim();
    % create the axes and set some properties
    set(ax,'box', 'off', 'nextplot', 'add', 'XMinorTick', 'off', 'YMinorTick', 'off' );
    
    ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/4;
    ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/4;

    % Set the tick direction
    ax.TickDir = 'out';
    % draw the plot to generate the undocumented vertex data var
    drawnow()
    % extract the x axis vertext data
    % X, Y and Z row of the start and end of the individual axle.
    vd = get(ax.XAxis.Axle,'VertexData');
    % reset the zero value
    vd(1,1) = xlims(1);
    % Update the vertex data
    set(ax.XAxis.Axle,'VertexData',vd);
    
    % repeat for Y (set 2nd row)
    vd = get(ax.YAxis.Axle,'VertexData');
    vd(2,1) = ylims(1);
    set(ax.YAxis.Axle,'VertexData',vd);
%     set(gca,'color','none')

    %     % extract the x axis vertext data
    %     % X, Y and Z row of the start and end of the individual axle.
    %     ax.XAxis.Axle.VertexData(1,1) = 0;
    %     % repeat for Y (set 2nd row)
    %     ax.YAxis.Axle.VertexData(2,1) = 0;
    %
    %     % You can modify the minor Tick values by modifying the vertex data
    %     % for them, e.g. remove any minor ticks below 0
    %     ax.XAxis.MinorTickChild.VertexData(:,ax.XAxis.MinorTickChild.VertexData(1,:)<0) = [];
    %     ax.YAxis.MinorTickChild.VertexData(:,ax.YAxis.MinorTickChild.VertexData(2,:)<0) = [];
end

end