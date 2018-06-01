function [lc, quietZoneIndex] = tb_investigateZone( raw, samples, fs, ecg_width, qz_width )
%               define some variables for the filter generator
                tb_lp = 30; % width of transition band [Hz]
                fc1_lp = 60; % location of cutoff freq [Hz]
                Rs_lp = 80; % stopband attentuaion [dB down]
                Rp_lp = 0.01; % passband attenuation ripple [dB]
                choice_lp = 2; % % 1 - butter, 2 - cheby1, 3 - cheby2, 4 - elliptic
                Imp_Length = 250; % filter gets cutoff after this many samples
                TOL_lp = 0.9;
                [ num, den ] = TB_LP_filterGenerator( fc1_lp, tb_lp, fs, Rs_lp, Rp_lp, Imp_Length, choice_lp);
                
                
                verify = filter(num, den, raw);

%               [pk,lc] = findpeaks(verify,samples,'SortStr','none', 'MinPeakHeight',TOL_lp*max(verify));
                [pk,lc] = findpeaks(verify,samples,'SortStr','none','MinPeakDistance', fs*5/8, 'MinPeakHeight',6);                
                if ( floor(lc(1) - ecg_width/2) < 0 )
                    pk = pk(2:end);
                    lc = lc(2:end);
                end
                
                if ( floor(lc(length(lc)) + ecg_width/2) > length(raw) )
                    pk = pk(1:(end-1));
                    lc = lc(1:(end-1));
                end
                
                
                
%                 [pk_raw,lc_raw] = findpeaks(raw,samples,'SortStr','none','MinPeakDistance', fs/2);
%                 figure
%                 plot (samples, raw);hold on;
%                 plot(samples, verify); 
% 
%                 plot(lc,pk,'x');
%                 for k = 1:1:length(pk)
%                     text(lc(k),pk(k), ['\leftarrow',num2str(lc(k)),'=n_{',num2str(k),'}'], 'Color','green')
%                 end
%                 plot(lc_raw,pk_raw,'x');
%                 for k = 1:1:length(pk_raw)
%                     text(lc_raw(k),pk_raw(k), ['\leftarrow',num2str(lc_raw(k)),'=n_{',num2str(k),'}'], 'Color','green')
%                 end
%                 hold off;
% 
%                 legend('Raw','Clarified');
%                 title('Clarified ECG Zones');

%                 disp('estimate ECG vector for data set');
%                 disp(lc);
                
                quietZoneIndex = zeros(1,length(lc)-1);
                
                for k = 1:1:(length(lc)-1)
                    quietZoneIndex(1,k) = (lc(k+1)+lc(k))/2;
                end
%                 disp('estimate QZ vector for data set');
%                 disp(quietZoneIndex);
end

