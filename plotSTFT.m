function [] = plotSTFT(STFT_M, STFT_fm, sig, sig_td_filt, sig_axis, fs, cutoff_frq , num, sig_type, yLo, yHi, printFileName, directory)

%     STFT_M = 512;
%     STFT_fm = 2048;
    
    selectMeshPlot = 0;
    
    [sigSTFT_T,sigSTFT_F,sigSTFT_SP] = STFT([sig_axis/fs; sig]', STFT_M, STFT_fm);
    [sigSTFT_td_T,sigSTFT_td_F,sigSTFT_td_SP] = STFT([sig_axis/fs; sig_td_filt]', STFT_M, STFT_fm);

    
    if(selectMeshPlot == 1)
    figure; 
    subplot(3,1,1);
    plot(sig_axis/fs,sig,'linewidth',2,'color',[0,0,0.2]);
    grid on; grid minor;
    b = legend('Raw');
    %title([printFileName,' Segmented ',sig_type,'-',num2str(num)]);
    xlabel('Time [s]');
    ylabel('Magnitude [mV]');
    ylim([-15 15]);
    grid on; grid minor;
    set(b,'color','k');
    ax = gca;
    tb_PlotStyle( ax );


    ax2 = subplot(3,1,2);   
    mesh(sigSTFT_T,sigSTFT_F,sigSTFT_SP); %,'EdgeColor','none');
    %title([printFileName,' STFT: ',sig_type,'-',num2str(num),' - Unfiltered in TD']);
    grid on;grid minor;
    xlabel('Time [s]');
    ylabel('Freq [Hz]');
    ylim([yLo yHi]);
    shading interp; 
    colorbar;
    ax = gca;
    tb_PlotStyle( ax );

    ax3 = subplot(3,1,3);    
    mesh(sigSTFT_td_T,sigSTFT_td_F,sigSTFT_td_SP); %,'EdgeColor','none');
    %title([printFileName, ' STFT: ',sig_type,'-',num2str(num),' After Highpass filter - Cutoff Freq: ',num2str(cutoff_frq),'[Hz]']);
    grid on;grid minor;
    xlabel('Time [s]');
    ylabel('Freq [Hz]');
    ylim([0 400]);
    colorbar;
    shading interp; 
    ax = gca;
    tb_PlotStyle( ax );

    colormap(ax2,jet)
    colormap(ax3,jet)
    
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     f = getframe(gcf);
%     imwrite(f.cdata, ['.\images_51617\',printFileName,sig_type,num2str(num),'_CF',num2str(cutoff_frq),'.png']); 

  %  tb_plotSaver( gcf, [printFileName,sig_type,num2str(num)], directory )
    
    else
    figure; 
    subplot(3,1,1);
    plot(sig_axis/fs,sig,'linewidth',2,'color',[0,0,0.2]);
    grid on; grid minor;
    b = legend('Raw');
    %title([printFileName,' Segmented ',sig_type,'-',num2str(num)]);
    xlabel('Time [s]');
    ylabel('Magnitude [mV]');
    ylim([-30 30]);
    grid on; grid minor;
    %set(b,'color','k');
    ax = gca;
    tb_PlotStyle( ax );


    ax2 = subplot(3,1,2);
    %contour(sigSTFT_T,sigSTFT_F,sigSTFT_SP);
    contourf(sigSTFT_T,sigSTFT_F,sigSTFT_SP,'EdgeColor','none');
    %title([printFileName,' STFT: ',sig_type,'-',num2str(num),' - Unfiltered in TD']);
    grid on;
    xlabel('Time [s]');
    ylabel('Freq [Hz]');
    ylim([0 400]);
    %shading interp; 
    colorbar;
    ax = gca;
    tb_PlotStyle( ax2 );

    ax3 = subplot(3,1,3); 
    %contour(sigSTFT_td_T,sigSTFT_td_F,sigSTFT_td_SP);
    contourf(sigSTFT_td_T,sigSTFT_td_F,sigSTFT_td_SP,'EdgeColor','none');
    %title([printFileName, ' STFT: ',sig_type,'-',num2str(num),' After Highpass filter - Cutoff Freq: ',num2str(cutoff_frq),'[Hz]']);
    grid on;
    xlabel('Time [s]');
    ylabel('Freq [Hz]');
    ylim([0 400]);
    colorbar;
    %shading interp; 
    ax = gca;
    tb_PlotStyle( ax3 );

    colormap(ax2,jet)
    colormap(ax3,jet)
    
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     f = getframe(gcf);
%     imwrite(f.cdata, ['.\images_51617\',printFileName,sig_type,num2str(num),'_CF',num2str(cutoff_frq),'.png']); 

    %tb_plotSaver( gcf, [printFileName,'-STFT-',sig_type,num2str(num)], directory )        
    
    end

end

