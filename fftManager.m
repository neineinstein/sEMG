function [ P1_out, P1_2_out, f ] = fftManager( Z, Z_td, fs, TOL, zoneName, printFileName, directory, loopIndex, fileID )
% need to get the input zone and plot each signal in a useful way
% will try a (1) subplot way and an (2) overplot way
plotOutputFlag = 1;
color1 = [128 0 0]/255;
color2 = [0   0 128]/255;
color3 = [190 190 190]/255;

%% (1) 
[m, n] = size(Z_td);
%zoneName = 'Test ';
% [m2, n2] = size(ZX);
f = fs*(0:floor((n/2)))/n; 
fft_out = zeros(m,n);
P1_out = zeros(m,floor((n-1)/2+1));

[o, p] = size(Z);
%zoneName = 'Test ';
% [m2, n2] = size(ZX);
% f = fs*(0:floor((p/2)))/p; 
fft_out2 = zeros(o,p);
P1_2_out = zeros(o,floor((p-1)/2+1));

%         figure
%     for j = 1:1:m
% %        subplot(ceil(sqrt(m)),ceil(sqrt(m)), j)
%        subplot(3,3, j) % this line can cause alot of errors...how to fix this???
%         fft_out(j,:) = fft(Z_td(j,:));
% 
%         P2 = abs(fft_out(j,:)/length(Z_td(j,:)));
% 
%         P1 = P2(1:length(Z_td(j,:))/2+1);
%         P1_out(j,:) = P2(1:length(Z_td(j,:))/2+1);
% 
%         P1(2:end-1) = 2*P1(2:end-1);
%         P1_out(j,2:end-1) = 2*P1_out(j,2:end-1);
% 
% 
%         [pk,lc] = findpeaks(P1_out(j,:),f,'SortStr','descend', 'MinPeakHeight',TOL*max(P1));
%         
%         % grab the 5 max peaks
%         if(length(pk) >= 5)
%             pk = pk(1:5);
%             lc = lc(1:5);
%         end
% 
%         plot(f,P1_out(j,:),'y');
%         hold on;
%         plot(lc,pk,'x');
% 
%         for k = 1:1:length(pk)
%             text(lc(k),pk(k), ['\leftarrow',num2str(lc(k)),'Hz'])
%         end
%         title([zoneName,':',num2str(j),':Single-Sided Amplitude Spectrum'])
%         xlabel('F [Hz]')
%         ylabel(['|Y_{',zoneName,'}(f)|'])
%         hold off;
%         plot_settings( [0,0,0] );
% 
%     end
%     tb_plotSaver( gcf, [printFileName,zoneName], directory )

    for j = 1:1:m

        fft_out(j,:) = fft(Z_td(j,:));
        fft_out2(j,:) = fft(Z(j,:));
        
        P2 = abs(fft_out(j,:)/length(Z_td(j,:)));
        P2_2 = abs(fft_out2(j,:)/length(Z(j,:)));
        
        P1 = P2(1:floor(length(Z_td(j,:))/2+1));
        P1_2 = P2_2(1:floor(length(Z(j,:))/2+1));
        
        P1_out(j,:) = P2(1:floor(length(Z_td(j,:))/2+1));
        P1_2_out(j,:) = P2_2(1:floor(length(Z(j,:))/2+1));


        P1(2:end-1) = 2*P1(2:end-1);
        P1_2(2:end-1) = 2*P1_2(2:end-1);
        
        P1_out(j,2:end-1) = 2*P1_out(j,2:end-1);
        P1_2_out(j,2:end-1) = 2*P1_2_out(j,2:end-1);
        
        [pk,lc] = findpeaks(P1_2_out(j,:),f,'SortStr','descend');
        if(length(pk) >= 5)
            pk = pk(1:5);
            lc = lc(1:5);
        end
        

        [pk2,lc2] = findpeaks(P1_out(j,:),f,'SortStr','descend');
        if(length(pk2) >= 5)
            pk2 = pk2(1:5);
            lc2 = lc2(1:5);
        end
        
        sort(lc);
        sort(lc2);
        
        
        
        if (plotOutputFlag == 1)
            
%         fprintf(fileID,[zoneName,'-',num2str(j),',  %3.0f, %3.0f, %3.0f, %3.0f, %3.0f \n'],lc(1),lc(2),lc(3),lc(4),lc(5));
%         fprintf(fileID,[zoneName,'-TD Filt-',num2str(j),', %3.0f, %3.0f, %3.0f, %3.0f, %3.0f \n'],lc2(1),lc2(2),lc2(3),lc2(4),lc2(5));
            
            
        figure
        subplot(2,1,1)
         %%begin unfilt
        plot(f,P1_2_out(j,:),'Color',color2);
        hold on;
        plot(lc,pk,'x');

        for k = 1:1:length(pk)
            text(lc(k),pk(k), ['\leftarrow',num2str(lc(k)),'Hz'])
        end
        %title([zoneName,':',num2str(j),':Single-Sided Amplitude Spectrum'])
        xlabel('F [Hz]')
        ylabel(['|Y_{',zoneName,'}(f)|'])
        legend('Unfiltered');
        hold off;
        plot_settings( [0,0,0] );
        %%end unfilt
       
        
        

      
        subplot(2,1,2)
 %% begin filtered 
        plot(f,P1_out(j,:),'Color',color2);
        hold on;
        plot(lc2,pk2,'x');

        for k = 1:1:length(pk2)
            text(lc2(k),pk2(k), ['\leftarrow',num2str(lc2(k)),'Hz'])
        end
       % title([zoneName,'-TD Filt:',num2str(j),':Single-Sided Amplitude Spectrum'])
        xlabel('F [Hz]')
        ylabel(['|Y_{',zoneName,'}(f)|'])
        legend('Filtered');
        hold off;
        plot_settings( [0,0,0] );
        
        
        %% end td filt
        
       
%         tb_plotSaver( gcf, [printFileName,'-FFT-',zoneName,num2str(j)], directory )
        end
    end

end

