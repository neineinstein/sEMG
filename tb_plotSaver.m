function [  ] = tb_plotSaver( gcf, title, directory )
%TB_PLOTSAVER save current plot to a directory with a title
%   steps

%   1.) meld directory/title as totalFilePath
    totalFilePath = fullfile(directory, strcat(title,'.png'));

%   2.) make the plot full screen
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    
%   3.) get the plot as an image
     f = getframe(gcf);
     
%   4.) write the image to a file
     imwrite(f.cdata, totalFilePath);

end

