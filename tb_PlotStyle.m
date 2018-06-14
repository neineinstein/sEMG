function [  ] = tb_PlotStyle( ax )
% takes ax input and sets plot styles
% ignore this line
    mycolor = [0 0 0];
   % b.Color = [0.8 0.8 0.8];
    ax.Color = mycolor;
    ax.GridColor = mycolor;
    ax.MinorGridColor = mycolor;
    ax.XColor  = mycolor;
    ax.YColor  = mycolor;
    ax.MinorGridLineStyle = ':'; 
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.5;
    ax.MinorGridAlpha = 0.5;
    set(ax,'fontweight','bold','FontSize',12,'linewidth',1)
    set(ax,'Color',[1 1 1])
end

