function  plot_settings( mycolor )

    ax = gca;
    grid on; %grid minor;
    ax.Color = mycolor;
    ax.GridColor = mycolor;
    ax.MinorGridColor = mycolor;
    ax.XColor  = mycolor;
    ax.YColor  = mycolor;
    ax.MinorGridLineStyle = ':'; 
    ax.GridLineStyle = ':';
    ax.GridAlpha = 0.5;
    ax.MinorGridAlpha = 0.5;
    set(ax,'fontweight','bold','FontSize',16,'linewidth',1)
    set(ax,'Color',[1 1 1])
    
    % minimize white space
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
end

