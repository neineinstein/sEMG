function tb_plotTotalFFT( f, z, plotName )
%tb_plotTotalFFT Plots all of the FFTs generated from a single patient
%   Detailed explanation goes here
 [~, ~, n] = size(z);
%  data = zeros(m, n);
   figure
%    title('TD Filtered QZ FFTs');
    for j=1:1:n
        data = z(:,:,j);
        data( ~any(data,2), : ) = [];  %rows
        [m, ~] = size(data);
        for k=1:1:m
            %plot(f,z(k,:,j));
            plot(f,data(k,:)/max(data(k,:)));
            hold on;
        end
    end
    hold off;
    title(plotName);
    xlabel('F [Hz]')
    ylabel(['Magnitude'])
    plot_settings( [0, 0, 0] );
end

