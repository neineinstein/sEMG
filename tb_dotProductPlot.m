function tb_dotProductPlot( f, x, y , plotName)
    
    figure
    plot(f,(x.*y));
    title(plotName);
    xlabel('F [Hz]')
    ylabel(['Magnitude'])
    %ylabel(['|Y_{',zoneName,'}(f)|'])
    plot_settings( [0, 0, 0] );
    
end

