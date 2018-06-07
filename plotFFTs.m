function  plotFFTs( f_qz, f_ecg, raw_ECG_FFTs, filtered_ECG_FFTs, raw_QZ_FFTs, filtered_QZ_FFTs, printFileName, directory )

    %% plot all of the FFT signals
%     plotName = 'filtered QZ FFTs'; 
%     tb_plotTotalFFT( f_qz , filtered_QZ_FFTs, plotName);
% 
%     plotName = 'raw QZ FFTs'; 
%     tb_plotTotalFFT( f_qz , raw_QZ_FFTs, plotName);
% 
%     plotName = 'filtered ECG FFTs'; 
%     tb_plotTotalFFT( f_ecg , filtered_ECG_FFTs, plotName);
% 
%     plotName = 'raw ECG FFTs'; 
%     tb_plotTotalFFT( f_ecg , raw_ECG_FFTs, plotName);
%     Ignore this line
    doAverage = 0;

%% generate dot products against each other 
    [ dp_Filt_QZ_FFT ] = tb_zoneDotProduct( filtered_QZ_FFTs, doAverage);
    [ dp_Raw_QZ_FFT ] = tb_zoneDotProduct( raw_QZ_FFTs, doAverage );
    [ dp_Filt_ECG_FFT ] = tb_zoneDotProduct( filtered_ECG_FFTs, doAverage );
    [ dp_Raw_ECG_FFT ] = tb_zoneDotProduct( raw_ECG_FFTs, doAverage );
    
    plotName = 'Time Domain Filtered QZ and ECG Dot Product Plot';
    tb_dotProductPlot( f_ecg, dp_Filt_QZ_FFT, dp_Filt_ECG_FFT, plotName );
    plotName = 'QZ and ECG Dot Product Plot';
    tb_dotProductPlot( f_ecg, dp_Raw_QZ_FFT, dp_Raw_ECG_FFT, plotName );
    plotName = 'TD Filt QZ and Unfiltered ECG Dot Product Plot';
    tb_dotProductPlot( f_ecg, dp_Filt_QZ_FFT, dp_Raw_ECG_FFT, plotName );
    
    
   
    
    [ avg_filtered_QZ_FFTs ] = tb_avgZone( filtered_QZ_FFTs );
    [ avg_Raw_QZ_FFT ] = tb_avgZone( raw_QZ_FFTs );
    [ avg_Filt_ECG_FFT ] = tb_avgZone( filtered_ECG_FFTs );
    [ avg_Raw_ECG_FFT ] = tb_avgZone( raw_ECG_FFTs );
    
    %% i turned this off because it's hardcoded to my computer
    %directory = 'C:\---\sEMG_DotProductPlots';

    plotName = 'ECG against ECG average(dotProduct)';
    tb_dotProductPlot( f_ecg, avg_Raw_ECG_FFT, raw_ECG_FFTs(:,:,1), plotName ); 
    %tb_plotSaver( gcf, [printFileName,'ECGvsAvgECG'], directory );
    
    % compare 
    plotName = [printFileName,'Averaged TD Filt QZ and Raw ECG Dot Product Plot'];
    tb_dotProductPlot(f_ecg, avg_Raw_ECG_FFT, avg_filtered_QZ_FFTs, plotName);
    %tb_plotSaver( gcf, [printFileName,'Avg_FiltQZvsRawECG'], directory );
    
    plotName = [printFileName,'AvgRaw QZ and AvgFilt TD ECG Dot Product Plot'];
    tb_dotProductPlot(f_ecg, avg_Filt_ECG_FFT, avg_Raw_QZ_FFT, plotName);
    %tb_plotSaver( gcf, [printFileName,'Avg_FiltECGvsRawQZ'], directory );
    
end

