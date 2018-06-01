%%%% Bradley University %%% Yufeng Lu
%% revised for sEMG project 
%%% 9/7/2016
clear all; clc; close all;
for j = 1:1:1
%% prepare the data 
filename = 'as6m.mat';
file = matfile(filename, 'Writable', true);

f_sampling = 2048; %%
n_samples = f_sampling*5;
xls_begin = 9; %% start at B9
xls_end = n_samples-1+xls_begin; %% select n_samples
string_xls_name = ['103-normal-M-6-6-13.xlsx']; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS
string_sheet = ['R1_4'];   % sheet namey
string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
raw =xlsread(string_xls_name,string_sheet,string_cells);
%% 
% figure;
% plot(raw);grid on; title(['recorded data'  string_xls_name]);
% xlabel('# sample ');
% ylabel('Amplitude [mV]');

samples = 1:1:length(raw);

%% initial segmentation of ECG and quiet zones.  this needs to be improved by use of the pan tompkins algorithm
% These are break points between quiet zones and ECG pulses
b1 = 897;
b2 = 1072;
b3 = 2505;
b4 = 2674;
b5 = 4170;
b6 = 4346;
b7 = 5891;
b8 = 6082;
b9 = 7625;
b10 = 7807;
b11 = 9252;
b12 = 9432;
b13 = length(raw);
% qz1 = [1 897];

%qz1x = (1 <= samples) & (samples <= b1);
qz1x = samples(1:b1);
qz1 = raw(1:b1)';
% ecg1 = [898 1072];
ecg1x = samples(b1:b2);
ecg1 = raw(b1:b2)';
% qz2 = [1073 2505];
qz2x = samples(b2:b3);
qz2 = raw(b2:b3)';
% ecg2 = [2506 2632];
ecg2x = samples(b3:b4);
ecg2 = raw(b3:b4)';
% qz3 = [2633 4170];
qz3x = samples(b4:b5);
qz3 = raw(b4:b5)';
% ecg3 = [4171 4346];
ecg3x = samples(b5:b6);
ecg3 = raw(b5:b6)';
% qz4 = [4347 5891];
qz4x = samples(b6:b7);
qz4 = raw(b6:b7)';
% ecg4 = [5892 6082];
ecg4x = samples(b7:b8);
ecg4 = raw(b7:b8)';
% qz5 = [6083 7625];
qz5x = samples(b8:b9);
qz5 = raw(b8:b9)';
% ecg5 = [7626 7807];
ecg5x = samples(b9:b10);
ecg5 = raw(b9:b10)';
% qz6 = [7808 9252];
qz6x = samples(b10:b11);
qz6 = raw(b10:b11)';
% ecg6 = [9253 9432];
ecg6x = samples(b11:b12);
ecg6 = raw(b11:b12)';
% qz7 = [9433 length(raw)];
qz7x = samples(b12:b13);
qz7 = raw(b12:b13)';

figure
subplot(2,1,1)
plot(raw);grid on; title(['recorded data'  string_xls_name]);
grid on; grid minor;
xlabel('# sample ');
ylabel('Amplitude [mV]');

subplot(2,1,2)
hold on
plot(samples(ecg1x),ecg1,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(qz1x),qz1,'linewidth',1,'color',[1,0.8,0])
legend('ECG Segments','Quiet Zone Segments');
plot(samples(ecg2x),ecg2,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(ecg3x),ecg3,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(ecg4x),ecg4,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(ecg5x),ecg5,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(ecg6x),ecg6,'linewidth',1,'color',[0,0.7,0.2])
plot(samples(qz2x),qz2,'linewidth',1,'color',[1,0.8,0])
plot(samples(qz3x),qz3,'linewidth',1,'color',[1,0.8,0])
plot(samples(qz4x),qz4,'linewidth',1,'color',[1,0.8,0])
plot(samples(qz5x),qz5,'linewidth',1,'color',[1,0.8,0])
plot(samples(qz6x),qz6,'linewidth',1,'color',[1,0.8,0])
plot(samples(qz7x),qz7,'linewidth',1,'color',[1,0.8,0])
hold off;
grid on; grid minor;
title('Segmented Signal');
xlabel('# sample ');
ylabel('Amplitude [mV]');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
f = getframe(gcf);
imwrite(f.cdata, ['.\images\signals.png']); 


for j = 1:1
figure
subplot(2,3,1);
plot(ecg1);
title('ECG Pulse 1');

subplot(2,3,2); 
plot(ecg2);
title('ECG Pulse 2');

subplot(2,3,3);
plot(ecg3);
title('ECG Pulse 3');

subplot(2,3,4); 
plot(ecg4);
title('ECG Pulse 4');

subplot(2,3,5); 
plot(ecg5);
title('ECG Pulse 5');

subplot(2,3,6); 
plot(ecg6);
title('ECG Pulse 6');

figure
subplot(2,3,1);
plot(qz1);
title('Quiet Zone 1');

subplot(2,3,2); 
plot(qz2);
title('Quiet Zone 2');

subplot(2,3,3);
plot(qz3);
title('Quiet Zone 3');

subplot(2,3,4); 
plot(qz4);
title('Quiet Zone 4');

subplot(2,3,5); 
plot(qz5);
title('Quiet Zone 5');

subplot(2,3,6); 
plot(qz6);
title('Quiet Zone 6');

subplot(2,3,6); 
plot(qz7);
title('Quiet Zone 7');

end

% QZ = struct{one,zeros(
% 
% e = struct('one',{zeros(1,length(tf/T)-2)},'two',{zeros(1,length(tf/T)-2)},...
%     'three',{zeros(1,length(tf/T)-2)});

last_peak_loc = 7500; %% first 3 echoes
qrs_offset = round(0.1*f_sampling);
%%%%%%%%% test 3  (it takes longer to run)
preset_len = last_peak_loc-qrs_offset;

NEcho = 3;

% why the 1000 and 8000 here.  What happens if I change these numbers??
t_start = 1000;
t_end = 8000 ;
p = 2;
q = 11;
STFT_M = 256;%256
STFT_fm = f_sampling;

%% windowing by observation  
NSamples = t_end - t_start; %% even number

if (mod(NSamples,2) == 1) 
    NSamples = NSamples - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%
temp = raw(t_start:t_end)';
%%%%%%%%%%%%%%%%%%%%%%%
filter_sig = temp;
%%%%%%%%%%%%%%%%
temp = filter_sig;

fs = 2048;
delta = 1/fs; 
Ndata = length(temp);
%temp = temp./max(abs(temp));
%temp = temp - mean(temp); %% DC removal
original_sig = temp;

%%%%%%%%%%%%%%%%% EMD processing
x=original_sig;
%x = x - mean(x); 
%x = x./max(x);  
original_sig = x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
NSamples =length(x);
t = 1:NSamples;
% number of extrema and zero-crossings in residual
ner = NSamples;
nzr = NSamples;
r = x;
imf = [];
k = 1;
% iterations counter for extraction of 1 mode
Init_Niter=0;
% total iterations counter
NbIt=0;
INTERP = 'spline';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum number of iterations
MAXITERATIONS=500;
%main loop : requires at least 3 extrema to proceed
flag_stop_EMD=0;
%%%%if  flag_stop_EMD=1 then stop the EMD. 'cause there is no enough #extreme for EMD
while ~flag_stop_EMD 
    [indmin,indmax,indzer] = extr(r); %% get the index of extremes and zero-crossings
	ner = length(indmin) + length(indmax); %% # of extremes
	flag_stop_EMD = ner <3;	    
	% current mode
	m = r;
	% mode at previous iteration
	mp = m;    
    FIXE_H = 3;   %% stop_count 
    stop_count = 0;
	%%%%%%%%%%%%%%%%%
     try
        %%% iterpolation the signal and find the envelope of #min, #max
   		 [envmin,envmax,moyenne,nem,nzm] = envelope(t,m,INTERP);
   		if (abs(nzm-nem)>1)
			stop_sift = 0;
			stop_count = 0;
		else
			stop_count = stop_count+1;
            stop_count;
			stop_sift = (stop_count == FIXE_H);
		end
     catch        
% 		moyenne = zeros(1,length(m));
		stop_sift = 1;
    end    
    %%%%%%%%%%%%%%%%%%%    
    stop_sift;
	stop_count = 0;	    
   	while ~stop_sift && Init_Niter<MAXITERATIONS
        		%sifting
        	m = m - moyenne;
	       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             try        %%% iterpolation the signal and find the envelope of #min, #max
   		           [envmin,envmax,moyenne,nem,nzm] = envelope(t,m,INTERP);
   		        if (abs(nzm-nem)>1)
			         stop_sift = 0;
			         stop_count = 0;
                else
			        stop_count = stop_count+1;
                    stop_count;
			        stop_sift = (stop_count == FIXE_H);
                end
             catch                
		        moyenne = zeros(1,length(m));
		        stop_sift = 1;
             end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	    	mp = m;
        	Init_Niter=Init_Niter+1;
    	    NbIt=NbIt+1;
 	end % sifting loop
  	imf(k,:) = m;
  	%Final_Niter(k) = Init_Niter;
  	k = k+1;
 	r = r - m;
  	Init_Niter=0;
end %main loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot the EMD in time domain
% [NDecompose, junk] = size(imf);
%  NDecompose=4;  %%% changed
% figure;
%  for ii = 2:(NDecompose+1)
%        subplot(NDecompose+1,1,ii);
%        plot(imf((ii-1),:));
%        %plot([zeros(1,499) imf((ii-1),:)]);
%        xlim([1 NSamples]);
%        ylim([-1.2 1.2]);       
%        grid on;
%        if(ii == NDecompose +1) 
%           xlabel('Time');
%        end
%        
%  end   
%    subplot(NDecompose+1,1,1);
%    plot(sig_orig);grid on;
%    xlim([1 NSamples]);
%        ylim([-1.2 1.2]);
%    title('EMD results');    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
[N_imf N_samples ] = size(imf);

 
   fs = 2048;
   t1 = t.* (1/fs);
%  for ii = 1:(N_imf )
%       % subplot(N_imf+1,1,ii);
%       figure;
%        plot(t1,imf((ii),:));
%        xlim([0 max(t1)]);
%         ylim([-15 15]);  
%         title(['imf #' num2str(ii)]);
%        grid on;
%  end  
%    figure;
%     
%    plot(t1,sig_orig);grid on;
%    xlim([0 max(t1)]);
%    title('original data');
%      ylim([-15 15]); 
%      grid on;
  % title('Empirical Mode Decomposition(EMD) results');     
   
   
   
   
   
%    %%%%%% replot  (imf3 and imf4)
%   figure;
%    fs = 2048;
%    t1 = t.* (1/fs);
%     
%    subplot(3,1,1);
%    plot(t1, sig_orig);
%    xlim([0 max(t1)]);
%    ylim([-15 15]);
%       grid on;
%    ecg = imf(3,:)+imf(4,:)+imf(5,:);
%    
%    subplot(3,1,2);
%    plot(t1, ecg);
%    xlim([0 max(t1)]);
%   ylim([-15 15]); 
%       grid on;
% %    subplot(4,1,3);
% %    plot(t1, imf(4,:));
% %    xlim([0 max(t1)]);
% %    ylim([-1.2 1.2]); 
%    
%    subplot(3,1,3);
%    plot(t1, sig_orig - ecg);
%    xlim([0 max(t1)]);
%    %ylim([-1.2 1.2]); 
%      
%    ylim([-15 15]);
%    grid on;
end      
   
   fs = 2048;
   t1 = t.* (1/fs);
%     figure;
%  for ii = 2:6
%       subplot(6,1,ii);
%       
%        plot(t1,imf((ii-1),:));
%        xlim([0 max(t1)]);
%         ylim([-15 15]);  
%        % title(['imf #' num2str(ii)]);
%        grid on;
%  end  
   
    file.t1 = t1;
    file.original_sig = original_sig;
    
%      subplot(6,1,1);
%      plot(t1,original_sig);grid on;
%      xlim([0 max(t1)]);   
%      ylim([-15 15]); 
%      grid on;
%      title(' orignal, imf1, imf2, imf3, imf4, imf5 (top to down)');
          
     
     %% here's where Tony is working
        % ECG comes from the heart, gaps
        % Signals strong at 30Hz are ECG
        % Signals Strong at 100 to 300Hz are the eMD
        
        % sum the IMF's, generating "reconstructed signal"
        % Reconstructed is just ECG from IMFs
      est_IMF_ECG_sig = 0*imf(1,:);
      for i = p:q
        est_IMF_ECG_sig = est_IMF_ECG_sig + imf(i,:);
      end
      
      % generated the filtered signal
      est_EMD_sig = original_sig-est_IMF_ECG_sig;
      file.reconstructed_sig = est_IMF_ECG_sig;
      file.filtered_sig = est_EMD_sig;
      
      
      % break up IMF_ECG_ 
  
% qz1 does not exist
% IMF_ECG_qz1x = t1(1:b1);
% IMF_ECG_qz1 = est_IMF_ECG_sig(1:b1)';
% ecg1 = [898 1072];
IMF_ECG_ecg1x = t1(1:b2-t_start);
IMF_EMD_ecg1 = est_EMD_sig(1:b2-t_start);
IMF_ECG_ecg1 = est_IMF_ECG_sig(1:b2-t_start);
% qz2 = [1073 2505];
IMF_ECG_qz2x = t1(b2-t_start:b3-t_start);
IMF_EMD_qz2 = est_EMD_sig(b2-t_start:b3-t_start);
IMF_ECG_qz2 = est_IMF_ECG_sig(b2-t_start:b3-t_start);
% ecg2 = [2506 2632];
IMF_ECG_ecg2x = t1(b3-t_start:b4-t_start);
IMF_EMD_ecg2 = est_EMD_sig(b3-t_start:b4-t_start);
IMF_ECG_ecg2 = est_IMF_ECG_sig(b3-t_start:b4-t_start);
% qz3 = [2633 4170];
IMF_ECG_qz3x = t1(b4-t_start:b5-t_start);
IMF_EMD_qz3 = est_EMD_sig(b4-t_start:b5-t_start);
IMF_ECG_qz3 = est_IMF_ECG_sig(b4-t_start:b5-t_start);
% ecg3 = [4171 4346];
IMF_ECG_ecg3x = t1(b5-t_start:b6-t_start);
IMF_EMD_ecg3 = est_EMD_sig(b5-t_start:b6-t_start);
IMF_ECG_ecg3 = est_IMF_ECG_sig(b5-t_start:b6-t_start);
% qz4 = [4347 5891];
IMF_ECG_qz4x = t1(b6-t_start:b7-t_start);
IMF_EMD_qz4 = est_EMD_sig(b6-t_start:b7-t_start);
IMF_ECG_qz4 = est_IMF_ECG_sig(b6-t_start:b7-t_start);
% ecg4 = [5892 6082];
IMF_ECG_ecg4x = t1(b7-t_start:b8-t_start);
IMF_EMD_ecg4 = est_EMD_sig(b7-t_start:b8-t_start);
IMF_ECG_ecg4 = est_IMF_ECG_sig(b7-t_start:b8-t_start);
% qz5 = [6083 7625];
IMF_ECG_qz5x = t1(b8-t_start:end);
IMF_EMD_qz5 = est_EMD_sig(b8-t_start:end);
IMF_ECG_qz5 = est_IMF_ECG_sig(b8-t_start:end);
% % ecg5 = [7626 7807];
% IMF_ECG_ecg5x = t1(b9:b10);
% IMF_ECG_ecg5 = est_IMF_ECG_sig(b9:b10)';
% % qz6 = [7808 9252];
% IMF_ECG_qz6x = t1(b10:b11);
% IMF_ECG_qz6 = est_IMF_ECG_sig(b10:b11)';
% % ecg6 = [9253 9432];
% IMF_ECG_ecg6x = t1(b11:b12);
% IMF_ECG_ecg6 = est_IMF_ECG_sig(b11:b12)';
% % qz7 = [9433 length(raw)];
% IMF_ECG_qz7x = t1(b12:b13);
% IMF_ECG_qz7 = est_IMF_ECG_sig(b12:b13)';
      
      
%% design a filter to apply in the time domain
%  we're going to try to design a high pass filter to let frequencies above
%  80Hz through
for fsc1 = 40:5:100
close all;
% High Pass Elliptical Filter
delta_p = 0.01;
delta_s = 0.001; 
Rs = -20*log10(delta_s)
Rp = -20*log10(1 - delta_p)
Rs = 80;
Rp = .01;
%fsc1 = 20;
fc1 = fsc1+3;  
fp1 = fc1*2/fs;
fs1 = fsc1*2/fs;
Wp = [fp1];
Ws = [fs1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% [N_ellip,Wn_ellip] = ellipord(Wp,Ws,Rp,Rs);
% [num_ellip,dem_ellip] = ellip(N_ellip,Rp,Rs,Wp);      
      
 hpFilt = designfilt('highpassfir','PassbandFrequency',fp1, ...
         'StopbandFrequency',fs1,'PassbandRipple',Rp, ...
         'StopbandAttenuation',Rs,'DesignMethod','kaiserwin');     
      
      
 % plot the design filter     
  fvtool(hpFilt);
    f = getframe(gcf);
    imwrite(f.cdata, ['.\images\filter_cutoffHz_',num2str(fsc1),'.png']);  
 % Generate the impulse response for time domain filtering
 [h_highPass,t] = impz(hpFilt);
 [hw_highPass,w] = freqz(hpFilt);
 
 
%   [h_highPass,t] = impz(num_ellip, dem_ellip);
% show the signals we're looking at      
       y = 50;
%       figure;
%       subplot(3,1,1);
%       plot(t1,original_sig);
%       title('Original Signal')
%       grid on;
%       ylim([-y y])
%       subplot(3,1,2);
%       plot(t1,est_IMF_ECG_sig);
%       title(['Reconstructed Signal(estimated ECG) (imf' num2str(p) '+...+imf' num2str(q) ')']);
%       grid on;
%       ylim([-y y])
%       subplot(3,1,3);
%       plot(t1,est_EMD_sig);
%       title('Filtered Signaln(Estimated EMD)');
%       grid on;
%       ylim([-y y])
%      
% now apply the filter designed above and plot those results
% going to apply the filter excesively and plot the results to compare
%     est_ECG_after_TD_filt = filter(h_highPass,1, est_IMF_ECG_sig);
%     orig_sig_after_TD_filt = filter(h_highPass, 1, original_sig);
%     est_EMD_after_TD_filt = filter(h_highPass, 1, est_EMD_sig);
    
   % STFT_fm = 32;
%     STFT_M = 256;

    yLo = 0;
    yHi = 100;
    %% apply the filter to ECG signals individually
    %ECG 1
        plotComparisons(h_highPass, ecg1, ecg1x, IMF_ECG_ecg1, IMF_ECG_ecg1x, IMF_EMD_ecg1, fs, fsc1 ,1, 'ECG', yLo, yHi)    
     %ECG2
        plotComparisons(h_highPass, ecg2, ecg2x, IMF_ECG_ecg2, IMF_ECG_ecg2x,IMF_EMD_ecg2, fs, fsc1 ,2, 'ECG', yLo, yHi)    
     %ECG3
        plotComparisons(h_highPass, ecg3, ecg3x, IMF_ECG_ecg3, IMF_ECG_ecg3x,IMF_EMD_ecg3, fs, fsc1 ,3, 'ECG', yLo, yHi)  
     %ECG4
        plotComparisons(h_highPass, ecg4, ecg4x, IMF_ECG_ecg4, IMF_ECG_ecg4x,IMF_EMD_ecg4, fs, fsc1 ,4, 'ECG', yLo, yHi)  
%      %ECG5
%         plotComparisons(h_highPass, ecg5, ecg5x, fs, fsc1 ,5, 'ECG')   
%      %ECG6
%         plotComparisons(h_highPass, ecg6, ecg6x, fs, fsc1 ,6, 'ECG') 



    yLo = 0;
    yHi = 400;
     %QZ1   
%         plotComparisons(h_highPass, qz1, qz1x, IMF_ECG_qz1, IMF_ECG_qz1x, fs, fsc1 ,1, 'Quiet Zone')
     %QZ2   
        plotComparisons(h_highPass, qz2, qz2x, IMF_ECG_qz2, IMF_ECG_qz2x, IMF_EMD_qz2, fs, fsc1 ,2, 'Quiet Zone', yLo, yHi) 
     %QZ3   
        plotComparisons(h_highPass, qz3, qz3x, IMF_ECG_qz3, IMF_ECG_qz3x, IMF_EMD_qz3, fs, fsc1 ,3, 'Quiet Zone', yLo, yHi) 
     %QZ4   
        plotComparisons(h_highPass, qz4, qz4x, IMF_ECG_qz4, IMF_ECG_qz4x, IMF_EMD_qz4, fs, fsc1 ,4, 'Quiet Zone', yLo, yHi) 
%      %QZ5   
%         plotComparisons(h_highPass, qz5, qz5x, fs, fsc1 ,5, 'Quiet Zone') 
%      %QZ6   
%         plotComparisons(h_highPass, qz6, qz6x, fs, fsc1 ,6, 'Quiet Zone') 
%      %QZ7   
%         plotComparisons(h_highPass, qz7, qz7x, fs, fsc1 ,7, 'Quiet Zone') 
     
  end       
        
   % est_EMD_after_TD_filt = filter(h,1, est_EMD_sig);
%       figure
%       subplot(3,1,1);
%       plot(t1, est_ECG_after_TD_filt); hold on;
%       plot(t1, est_IMF_ECG_sig);
%       title(['ECG Signal After Highpass filter - Cutoff Freq: ',num2str(fsc1),'[Hz]']);
%       grid on;
%       ylim([-y y])      
%  
%       subplot(3,1,2);
%       plot(t1, orig_sig_after_TD_filt); hold on;
%       plot(t1, original_sig);
%       title(['Original Signal Highpass filter - Cutoff Freq: ',num2str(fsc1),'[Hz]']);
%       grid on;
%       ylim([-y y]) 
% 
%       subplot(3,1,3);
%       plot(t1, est_EMD_sig); hold on;
%       plot(t1, est_EMD_after_TD_filt);
%       title(['EMD Signal Highpass filter - Cutoff Freq: ',num2str(fsc1),'[Hz]']);
%       grid on;
%       ylim([-y y])       
%       
%       STFT_M = 512;
%       %% perform STFT on original signal, ECG estimate, and EMD estimate
%     [reconstructed_T,reconstructed_F,reconstructed_SP] = STFT([t1; est_IMF_ECG_sig]', STFT_M, STFT_fm);
%     [ECG_TD_FILT_T,ECG_TD_FILT_F,ECG_TD_FILT_SP] = STFT([t1; est_ECG_after_TD_filt]', STFT_M, STFT_fm);
%     
%     [filtered_T,filtered_F,filtered_SP] = STFT([t1; est_EMD_sig]', STFT_M, STFT_fm);
%     [est_EMD_after_TD_filt_T,est_EMD_after_TD_filt_F,est_EMD_after_TD_filt_SP] = STFT([t1; est_EMD_after_TD_filt]', STFT_M, STFT_fm);
%      
%     [original_T,original_F,original_SP] = STFT([t1; original_sig]', STFT_M, STFT_fm);
%     [orig_sig_after_TD_filt_T,orig_sig_after_TD_filt_F,orig_sig_after_TD_filt_SP] = STFT([t1; orig_sig_after_TD_filt]', STFT_M, STFT_fm);
 
    
    
    %% PLOT COMPARISONS
    
    % subplot to show ECG1 as
    % 1 segmented
    % 2 clipped version from total filtered
    % 3 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     reconstructed_SP = filter(hw_highPass, 1, reconstructed_SP);
%     % verify that the filter works by plotting STFT of TD filtered and the
%     % FD filtered STFT of a signal
%     % compare the EMD vs the filtered here
%     figure; 
%     subplot(2,1,1);
%     contour(ECG_TD_FILT_T,ECG_TD_FILT_F,ECG_TD_FILT_SP);
%     title(sprintf('TD Filtered then STFTed Estimated ECG Signal'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 400]);
% 
%     subplot(2,1,2);
%     contour(reconstructed_T,reconstructed_F,reconstructed_SP*abs(hw_highPass));
%     title(sprintf('STFT ECG Estimate then Filter in FD'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 400]);   
    
    
    
    
%     % compare the EMD vs the filtered here
%     figure; 
%     subplot(2,1,1);
%     contour(filtered_T,filtered_F,filtered_SP);
%     title(sprintf('EMD Signal From IMF Subtraction'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 400]);
% 
%     subplot(2,1,2);
%     contour(ECG_TD_FILT_T,ECG_TD_FILT_F, ECG_TD_FILT_SP);
%     title(sprintf('STFT of TD Filtered IMF ECG Estimate'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 400]);    
%     
%       % compare the ECG filtering here
%     figure; 
%     subplot(2,1,1);
%     contour(reconstructed_T,reconstructed_F,reconstructed_SP);
%     title(sprintf('STFT of IMF ECG Estimated'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 100]);
% 
%     subplot(2,1,2);
%     contour(ECG_TD_FILT_T,ECG_TD_FILT_F, ECG_TD_FILT_SP);
%     title(sprintf('STFT of TD Filtered IMF ECG Estimate'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 400]);
%     
    
    
    
    
    
    
%    % heres the segment to STFT
%       tstart = 0.75;
%       tstop = 0.85;
%       xVector = t1(floor(tstart*fs):floor(tstop*fs));  % x axis data 
%       yVector = original_sig(floor(tstart*fs):floor(tstop*fs)); % y axis data
%       yVector_rec = est_IMF_ECG_sig(floor(tstart*fs):floor(tstop*fs)); % y axis reconstructed 
%       yVector_filt = est_EMD_sig(floor(tstart*fs):floor(tstop*fs)); % y axis filtered 
%       seg_original = [xVector;yVector];
%       seg_rec = [xVector;yVector_rec];
%       seg_filt = [xVector;yVector_filt];
%       figure
%       plot( seg_original(1,:), seg_original(2,:));
%       title('ECG Segment original')
%       grid on; grid minor;
%       xlabel('time [s]');
%       xlim([0.75 0.85])
%       ylabel('Volts [mV]');
%       figure
%       plot( seg_rec(1,:), seg_rec(2,:));
%       title('ECG Segment reconstructed')
%       grid on; grid minor;
%       xlabel('time [s]');
%       xlim([0.75 0.85])
%       ylabel('Volts [mV]');
%       figure
%       plot( seg_filt(1,:), seg_filt(2,:));
%       title('ECG Segment filtered')
%       grid on; grid minor;
%       xlabel('time [s]');
%       xlim([0.75 0.85])
%       ylabel('Volts [mV]');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STFT_F = 64;   
%     [seg_original_T,seg_original_F,seg_original_SP] = STFT([seg_original(1,:); seg_original(2,:)]', STFT_F, STFT_fm);
%     [seg_rec_T,seg_rec_F,seg_rec_SP] = STFT([seg_rec(1,:); seg_rec(2,:)]', STFT_F, STFT_fm);
%     [seg_filt_T,seg_filt_F,seg_filt_SP] = STFT([seg_filt(1,:); seg_filt(2,:)]', STFT_F, STFT_fm);
%     figure; 
%     subplot(3,1,1);
%     contour(seg_original_T,seg_original_F,seg_original_SP);
%     title(sprintf('STFT of original signal'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 200]);
%    
%     subplot(3,1,2);
%     contour(seg_rec_T,seg_rec_F, seg_rec_SP);
%     title(sprintf('STFT of reconstructed signal'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 200]);
%     
%     subplot(3,1,3);
%     contour(seg_filt_T,seg_filt_F, seg_filt_SP);
%     title(sprintf('STFT of filtered signal'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 200]);
%     grid on;
      
      
% % %% Design the filter 
% %     
% % Fs = 2048;  % Sampling Frequency
% % Fstop = 100;      % Stopband Frequency
% % Fpass = 300;     % Passband Frequency
% % Astop = 80;      % Stopband Attenuation (dB)
% % Apass = 1;       % Passband Ripple (dB)
% % match = 'both';  % Band to match exactly
% %  
% % % Construct an FDESIGN object and call its ELLIP method.
% % h  = fdesign.bandpass(Fstop, Fpass, Astop, Apass, Fs);
% % Hd = design(h, 'ellip', 'MatchExactly', match);
%  
% %Apply filter design to ECG Segment
% y_sig = est_ECG_after_TD_filt;
% 
% [y_sig_T,y_sig_F,y_sig_SP] = STFT([t1; y_sig]', STFT_F, STFT_fm);
% 
%  figure; 
%     
%     contour(y_sig_T,y_sig_F,y_sig_SP);
%     title(sprintf('FD STFT of original signal'));
%     grid on;grid minor;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 200]);
%  
% y_sig_filtered = filter(h,sig_sec_filtered);
% [y_sig_filtered_T,y_sig_filtered_F,y_sig_filtered_SP] = STFT([t3; y_sig_filtered]', STFT_F, STFT_fm);
%  
% y_sig_reconstructed = filter(h,sig_sec_imf);
% [y_sig_reconstructed_T,y_sig_reconstructed_F,y_sig_reconstructed_SP] = STFT([t3;y_sig_reconstructed]', STFT_F, STFT_fm);
%  
% % Apply filter design to the entire signal
% entire_y = filter (h,original_sig);
% [entire_y_T,entire_y_F,entire_y_SP] = STFT ([t1; entire_y]',STFT_M,STFT_fm);
%  
% entire_y_reconstructed = filter (h,est_IMF_ECG_sig);
% [entire_y_reconstructed_T,entire_y_reconstructed_F,entire_y_reconstructed_SP] = STFT ([t1; entire_y_reconstructed]',STFT_M,STFT_fm);
%  
% entire_y_filtered = filter (h,est_EMD_sig);
% [entire_y_filtered_T,entire_y_filtered_F,entire_y_filtered_SP] = STFT ([t1; entire_y_filtered]',STFT_M,STFT_fm);
% 
% %% plot the signals with filter applied
% for j = 1:1:1
% figure; 
%     subplot(3,1,1);
%     contour(A_original_T,A_original_F,A_original_SP);
%     title(sprintf('ellip Filter original signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     
%     subplot(3,1,2);
%     contour(A_reconstructed_sig_T,A_reconstructed_sig_F,A_reconstructed_sig_SP);
%     title(sprintf('ellip Filter reconstructed signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     
%     subplot(3,1,3);
%     contour(A_filtered_sig_T,A_filtered_sig_F,A_filtered_sig_SP);
%     title(sprintf('ellip Filter filtered signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     end  
% %% Plot Contour    
% for j = 1:1:1
% figure; 
%     subplot(3,1,1);
%     contour(entire_y_T,entire_y_F,entire_y_SP);
%     title(sprintf('Design Filter original signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     
%     subplot(3,1,2);
%     contour(entire_y_reconstructed_T,entire_y_reconstructed_F,entire_y_reconstructed_SP);
%     title(sprintf('Design Filter reconstructed signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     
%     subplot(3,1,3);
%     contour(entire_y_filtered_T,entire_y_filtered_F,entire_y_filtered_SP);
%     title(sprintf('Design Filter filtered signal'));
%     grid on;
%     xlabel('Time [S]');
%     ylabel('Freq [Hz]');
%     ylim([0 500]);
%     grid on;
%     
%     % Plot in time domain
%     figure;
%     subplot(3,1,1);
%     plot( A_original_sig);
%     title('Design filter Original signal in time domain')
%     grid on
%     
%     subplot(3,1,2);
%     plot(entire_y_reconstructed);
%     title('Design filter reconstructed signal in time domain')
%     grid on
%    
%     subplot(3,1,3);
%     plot( entire_y_filtered);
%     title('Design filter filtered signal in time domain')
%     grid on
%     
% end