%figure;plot(abs(fft(sig_residue)));grid on


function [outputArg1,outputArg2] = residue_difference(pat1, pat2)
% Author: Alex Topping
%   Function shows the average and difference of two (todo: or more)
%   patient data files
%   Where x is the name patient 1 xcel data
%   Where y is the name of patient 2 xcel data


%%%% Bradley University %%% Yufeng Lu
%% revised for sEMG project 
%%% created: 9/7/2016
%%% revised:  5/24/2018 



%% prepare the data 

%%% initial:  when the subject came to the lab
%%% 10-min resting time: 10 mins after the initial capture
%% for initial or 10-min resting,  20 measurements are done. 10 for right side, 10 for left side
%% For example: R1_5:   right side, initial data capture,  the 5th data set


 filename = 'init_test.mat';  %% need change the file name so that all results can be saved.
 
 file = matfile(filename, 'Writable', true);


%%%%%%%%% 
f_sampling = 2048; %%sampling rate
n_samples = f_sampling*5; %% pick 5 second long data
xls_begin = 9; %% start at sheet B9
xls_end = n_samples-1+xls_begin; %% select n_samples

%string_xls_name = ['103-normal-M-6-6-13.xlsx']; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS

string_xls_name = [pat1]; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS


%string_sheet = ['R1_5'];   % get the data from R1_5 sheet  

string_sheet = ['R1_6'];   % get the data from R1_6 sheet  

%string_sheet = ['R1_7'];   % get the data from R1_6 sheet  


string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
raw = -xlsread(string_xls_name,string_sheet,string_cells); %% 

figure;
plot(raw);grid on; title(['recorded data'  string_xls_name]);
xlabel('# sample ');
ylabel('Amplitude [mV]');



%%% by observation
%%% choose t_start to avoid starting with incomplete dominant echo
%%% choose t_end to avoid starting with incomplete dominant echo 
%%% make sure there is no obvious active muscle signal. (greater than 5 - 10 mv)

%%%% Set a breakpoint here to manually set t_start, t_end, and NEcho

%%%R1-5  103-normal-M-6-6-13
% t_start = 800;
% t_end = 7500;
% NEcho = 4; %% not used in the EMD process  (by observation)


% %%  R1-6   103-normal-M-6-6-13
% t_start =1;
% t_end = 10000;
% NEcho = 6; %% not used in the EMD process  (by observation)


% %%  R1-7   103-normal-M-6-6-13
% t_start =1;
% t_end = 6000;
% NEcho = 4; %% not used in the EMD process  (by observation)
% 


%%  R1-6   (203F-AS)
t_start =400;
t_end = 10000;
NEcho = 4; %% not used in the EMD process  (by observation)


%%% 




%% set the NSamples  
NSamples = t_end - t_start; %% even number
if (mod(NSamples,2) == 1) 
    NSamples = NSamples - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%
temp = raw(t_start:t_end)';
% temp = temp./max(abs(temp));
temp = temp - mean(temp); %% DC removal
%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% EMD processing

NSamples =length(temp);

if (mod(NSamples,2) == 1) 
    NSamples = NSamples - 1;
end

x=temp(1:NSamples);



original_sig = x;

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
NbIt=0; %% counter to track the # of iteration 
INTERP = 'spline';%% interpolation function used 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAXITERATIONS=500; % maximum number of iterations
%main loop : requires at least 3 extrema to proceed
flag_stop_EMD=0;% flag used to stop EMD process. (no enough #extreme for interpolation)

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
        %%% interpolate the signal and find the envelope of #min, #max
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
   	while ~stop_sift & Init_Niter<MAXITERATIONS
        		%sifting
        	m = m - moyenne;
	       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             try        %%% interpolate the signal and find the envelope of #min, #max
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
[NDecompose junk] = size(imf);
    
   fs = 2048;
   t1 = t.* (1/fs);
    figure;
    
    %% plot IMF 1 through 5 
 for ii = 2:6
      subplot(6,1,ii);
      
       plot(t1,imf((ii-1),:));
       xlim([0 max(t1)]);
        ylim([-15 15]);  
       % title(['imf #' num2str(ii)]);
       grid on;
 end  
   
    file.t1 = t1;
    file.original_sig = original_sig;
    
%%5 plot the original signal
    subplot(6,1,1);
   plot(t1,original_sig);grid on;
   xlim([0 max(t1)]);   
     ylim([-15 15]); 
     grid on;
     title(' orignal, imf1, imf2, imf3, imf4, imf5 (top to down)');
     
     %%%%%%%%%%%%%%%%%%%%%
         figure;
		 title(' imf6 to 11 (top to down)');
		 
         
 %%% plot IMF 6 through 11
 for ii = 1:5
      subplot(6,1,ii);
      
       plot(t1,imf((ii+5),:));
       xlim([0 max(t1)]);
        ylim([-15 15]);  
       % title(['imf #' num2str(ii)]);
       grid on;
 end  
      
     
 %%%%%%%
 
 eng_imf = zeros(NDecompose, 1);
 max_imf =  zeros(NDecompose, 1);
 
  var_imf =  zeros(NDecompose, 1);
 
 for ii = 1 : NDecompose
    
     
 eng_imf(ii) = sum(imf(ii,:).^2);    
     
 max_imf(ii) = max(abs(imf(ii,:)));    
     
 var_imf(ii) = var(imf(ii,:));  
 end
 
 
	 
	 %%%
	 p = 3;  % pick IMF 3 through 10
	 q = 5; 
	 
      reconstructed_sig_EMD = 0*imf(1,:);
      for i = p:q
        reconstructed_sig_EMD = reconstructed_sig_EMD + imf(i,:);
      end
	  
      filtered_sig = original_sig-reconstructed_sig_EMD;      
      file.reconstructed_sig = reconstructed_sig_EMD;
      file.filtered_sig = filtered_sig;
      
      y = 50;
      figure;
      subplot(3,1,1);
      plot(t1,original_sig);
      title('Original Signal')
      grid on;
      ylim([-y y])
      subplot(3,1,2);
      plot(t1,reconstructed_sig_EMD);
      title(['Reconstructed Signal (imf' num2str(p) '+...+imf' num2str(q) ')']);
      grid on;
      ylim([-y y])
      subplot(3,1,3);
      plot(t1,filtered_sig);
      title('Filtered Signal');
      grid on;
      ylim([-y y])
     

    
    

   
 
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  reconstructed_sig_EMD
    
    
sig_seg = zeros(NEcho,NSamples);
window_size = round(0.1 *f_sampling); %%% 0.1 second peak offset 

sig_residue = reconstructed_sig_EMD; %%%%%%%%%%%%%%
sig  = reconstructed_sig_EMD; %%%%%%%%%%%%%%

peak_loc_seg = zeros(NEcho,1);


figure;
plot(reconstructed_sig_EMD);
hold on;
%%% segment the data
for ii = 1 : NEcho
   [peak_val  peak_loc] = max(abs(sig_residue));   
   peak_loc_seg(ii) = peak_loc;
   peak_begin = max(1,peak_loc - window_size);
   peak_end = min(NSamples, peak_loc + window_size);
   sig_seg(ii,peak_begin:peak_end) = sig_residue(peak_begin:peak_end);    
   sig_residue(peak_begin:peak_end)= 0;
   plot(sig_seg(ii,:),'r');
end
title('segementated signal'); grid on;






figure;subplot(211); 
plot(original_sig);grid on; hold on; plot(reconstructed_sig_EMD,'g');grid on; plot(sig_residue,'r');
subplot(212);
plot(abs(fft(original_sig))); hold on; plot(abs(fft(reconstructed_sig_EMD)),'g'); grid on; hold on; plot(abs(fft(sig_residue)),'r');
xlim([0 2000]);



figure;subplot(311); 
plot(original_sig);grid on; subplot(312); plot(reconstructed_sig_EMD,'g');grid on; subplot(313); plot(sig_residue,'r'); grid on;

figure;
subplot(311);
plot(abs(fft(original_sig))); grid on;  xlim([0 2000]);  subplot(312); plot(abs(fft(reconstructed_sig_EMD)),'g'); grid on; xlim([0 2000]); subplot(313); plot(abs(fft(sig_residue)),'r');
xlim([0 2000]); grid on;



STFT_M = 256;
STFT_fm = f_sampling;


    [orig_T,orig_F,orig_SP] = STFT([t1; original_sig]', STFT_M, STFT_fm);
    figure; 
    contour(orig_T,orig_F,orig_SP);
    title(sprintf('Contour of STFT of Orignal Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
    
    file.orig_T = orig_T;
    file.orig_F = orig_F;
    file.orig_SP = orig_SP;


    [reconstructed_T,reconstructed_F,reconstructed_SP] = STFT([t1; reconstructed_sig_EMD]', STFT_M, STFT_fm);
    figure; 
    contour(reconstructed_T,reconstructed_F,reconstructed_SP);
    title(sprintf('Contour of STFT of EMD Reconstructed Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
    
    file.reconstructed_T = reconstructed_T;
    file.reconstructed_F = reconstructed_F;
    file.reconstructed_SP = reconstructed_SP;
    
    [residue_T,residue_F, residue_SP] = STFT([t1; sig_residue]', STFT_M, STFT_fm);
    figure; 
    contour(residue_T,residue_F,residue_SP);
    title(sprintf('Contour of STFT of Residue Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
  
    file.filtered_T = residue_T;
    file.filtered_F = residue_F;
    file.filtered_SP = residue_SP;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    removed_sig = original_sig - reconstructed_sig_EMD; 
    
   [removed_sig_T,removed_sig_F, removed_sig_SP] = STFT([t1; removed_sig]', STFT_M, STFT_fm);
    figure; 
    contour(removed_sig_T,removed_sig_F, removed_sig_SP);
    title(sprintf('Contour of STFT of Residue Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
  
    file.filtered_T = removed_sig_T;
    file.filtered_F = removed_sig_F;
    file.filtered_SP = removed_sig_SP;
    
close all
pat_1_residue = abs(fft(sig_residue));
figure;plot(pat_1_residue);grid on
title("yo")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%% 
f_sampling = 2048; %%sampling rate
n_samples = f_sampling*5; %% pick 5 second long data
xls_begin = 9; %% start at sheet B9
xls_end = n_samples-1+xls_begin; %% select n_samples

%string_xls_name = ['103-normal-M-6-6-13.xlsx']; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS

string_xls_name = [pat2]; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS


%string_sheet = ['R1_5'];   % get the data from R1_5 sheet  

string_sheet = ['R1_6'];   % get the data from R1_6 sheet  

%string_sheet = ['R1_7'];   % get the data from R1_6 sheet  


string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
raw = -xlsread(string_xls_name,string_sheet,string_cells); %% 

figure;
plot(raw);grid on; title(['recorded data'  string_xls_name]);
xlabel('# sample ');
ylabel('Amplitude [mV]');



%%% by observation
%%% choose t_start to avoid starting with incomplete dominant echo
%%% choose t_end to avoid starting with incomplete dominant echo 
%%% make sure there is no obvious active muscle signal. (greater than 5 - 10 mv)

%%%% Set a breakpoint here to manually set t_start, t_end, and NEcho

%%%R1-5  103-normal-M-6-6-13
% t_start = 800;
% t_end = 7500;
% NEcho = 4; %% not used in the EMD process  (by observation)


% %%  R1-6   103-normal-M-6-6-13
% t_start =1;
% t_end = 10000;
% NEcho = 6; %% not used in the EMD process  (by observation)


% %%  R1-7   103-normal-M-6-6-13
% t_start =1;
% t_end = 6000;
% NEcho = 4; %% not used in the EMD process  (by observation)
% 


%%  R1-6   (203F-AS)
t_start =400;
t_end = 10000;
NEcho = 4; %% not used in the EMD process  (by observation)


%%% 




%% set the NSamples  
NSamples = t_end - t_start; %% even number
if (mod(NSamples,2) == 1) 
    NSamples = NSamples - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%
temp = raw(t_start:t_end)';
% temp = temp./max(abs(temp));
temp = temp - mean(temp); %% DC removal
%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%% EMD processing

NSamples =length(temp);

if (mod(NSamples,2) == 1) 
    NSamples = NSamples - 1;
end

x=temp(1:NSamples);



original_sig = x;

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
NbIt=0; %% counter to track the # of iteration 
INTERP = 'spline';%% interpolation function used 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAXITERATIONS=500; % maximum number of iterations
%main loop : requires at least 3 extrema to proceed
flag_stop_EMD=0;% flag used to stop EMD process. (no enough #extreme for interpolation)

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
        %%% interpolate the signal and find the envelope of #min, #max
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
   	while ~stop_sift & Init_Niter<MAXITERATIONS
        		%sifting
        	m = m - moyenne;
	       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
             try        %%% interpolate the signal and find the envelope of #min, #max
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
[NDecompose junk] = size(imf);
    
   fs = 2048;
   t1 = t.* (1/fs);
    figure;
    
    %% plot IMF 1 through 5 
 for ii = 2:6
      subplot(6,1,ii);
      
       plot(t1,imf((ii-1),:));
       xlim([0 max(t1)]);
        ylim([-15 15]);  
       % title(['imf #' num2str(ii)]);
       grid on;
 end  
   
    file.t1 = t1;
    file.original_sig = original_sig;
    
%%5 plot the original signal
    subplot(6,1,1);
   plot(t1,original_sig);grid on;
   xlim([0 max(t1)]);   
     ylim([-15 15]); 
     grid on;
     title(' orignal, imf1, imf2, imf3, imf4, imf5 (top to down)');
     
     %%%%%%%%%%%%%%%%%%%%%
         figure;
		 title(' imf6 to 11 (top to down)');
		 
         
 %%% plot IMF 6 through 11
 for ii = 1:5
      subplot(6,1,ii);
      
       plot(t1,imf((ii+5),:));
       xlim([0 max(t1)]);
        ylim([-15 15]);  
       % title(['imf #' num2str(ii)]);
       grid on;
 end  
      
     
 %%%%%%%
 
 eng_imf = zeros(NDecompose, 1);
 max_imf =  zeros(NDecompose, 1);
 
  var_imf =  zeros(NDecompose, 1);
 
 for ii = 1 : NDecompose
    
     
 eng_imf(ii) = sum(imf(ii,:).^2);    
     
 max_imf(ii) = max(abs(imf(ii,:)));    
     
 var_imf(ii) = var(imf(ii,:));  
 end
 
 
	 
	 %%%
	 p = 3;  % pick IMF 3 through 10
	 q = 5; 
	 
      reconstructed_sig_EMD = 0*imf(1,:);
      for i = p:q
        reconstructed_sig_EMD = reconstructed_sig_EMD + imf(i,:);
      end
	  
      filtered_sig = original_sig-reconstructed_sig_EMD;      
      file.reconstructed_sig = reconstructed_sig_EMD;
      file.filtered_sig = filtered_sig;
      
      y = 50;
      figure;
      subplot(3,1,1);
      plot(t1,original_sig);
      title('Original Signal')
      grid on;
      ylim([-y y])
      subplot(3,1,2);
      plot(t1,reconstructed_sig_EMD);
      title(['Reconstructed Signal (imf' num2str(p) '+...+imf' num2str(q) ')']);
      grid on;
      ylim([-y y])
      subplot(3,1,3);
      plot(t1,filtered_sig);
      title('Filtered Signal');
      grid on;
      ylim([-y y])
     

    
    

   
 
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  reconstructed_sig_EMD
    
    
sig_seg = zeros(NEcho,NSamples);
window_size = round(0.1 *f_sampling); %%% 0.1 second peak offset 

sig_residue = reconstructed_sig_EMD; %%%%%%%%%%%%%%
sig  = reconstructed_sig_EMD; %%%%%%%%%%%%%%

peak_loc_seg = zeros(NEcho,1);


figure;
plot(reconstructed_sig_EMD);
hold on;
%%% segment the data
for ii = 1 : NEcho
   [peak_val  peak_loc] = max(abs(sig_residue));   
   peak_loc_seg(ii) = peak_loc;
   peak_begin = max(1,peak_loc - window_size);
   peak_end = min(NSamples, peak_loc + window_size);
   sig_seg(ii,peak_begin:peak_end) = sig_residue(peak_begin:peak_end);    
   sig_residue(peak_begin:peak_end)= 0;
   plot(sig_seg(ii,:),'r');
end
title('segementated signal'); grid on;






figure;subplot(211); 
plot(original_sig);grid on; hold on; plot(reconstructed_sig_EMD,'g');grid on; plot(sig_residue,'r');
subplot(212);
plot(abs(fft(original_sig))); hold on; plot(abs(fft(reconstructed_sig_EMD)),'g'); grid on; hold on; plot(abs(fft(sig_residue)),'r');
xlim([0 2000]);



figure;subplot(311); 
plot(original_sig);grid on; subplot(312); plot(reconstructed_sig_EMD,'g');grid on; subplot(313); plot(sig_residue,'r'); grid on;

figure;
subplot(311);
plot(abs(fft(original_sig))); grid on;  xlim([0 2000]);  subplot(312); plot(abs(fft(reconstructed_sig_EMD)),'g'); grid on; xlim([0 2000]); subplot(313); plot(abs(fft(sig_residue)),'r');
xlim([0 2000]); grid on;



STFT_M = 256;
STFT_fm = f_sampling;


    [orig_T,orig_F,orig_SP] = STFT([t1; original_sig]', STFT_M, STFT_fm);
    figure; 
    contour(orig_T,orig_F,orig_SP);
    title(sprintf('Contour of STFT of Orignal Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
    
    file.orig_T = orig_T;
    file.orig_F = orig_F;
    file.orig_SP = orig_SP;


    [reconstructed_T,reconstructed_F,reconstructed_SP] = STFT([t1; reconstructed_sig_EMD]', STFT_M, STFT_fm);
    figure; 
    contour(reconstructed_T,reconstructed_F,reconstructed_SP);
    title(sprintf('Contour of STFT of EMD Reconstructed Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
    
    file.reconstructed_T = reconstructed_T;
    file.reconstructed_F = reconstructed_F;
    file.reconstructed_SP = reconstructed_SP;
    
    [residue_T,residue_F, residue_SP] = STFT([t1; sig_residue]', STFT_M, STFT_fm);
    figure; 
    contour(residue_T,residue_F,residue_SP);
    title(sprintf('Contour of STFT of Residue Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
  
    file.filtered_T = residue_T;
    file.filtered_F = residue_F;
    file.filtered_SP = residue_SP;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    removed_sig = original_sig - reconstructed_sig_EMD; 
    
   [removed_sig_T,removed_sig_F, removed_sig_SP] = STFT([t1; removed_sig]', STFT_M, STFT_fm);
    figure; 
    contour(removed_sig_T,removed_sig_F, removed_sig_SP);
    title(sprintf('Contour of STFT of Residue Signal'));
    grid on;
    xlabel('Time [S]');
    ylabel('Freq [Hz]');
    grid on;
  
    file.filtered_T = removed_sig_T;
    file.filtered_F = removed_sig_F;
    file.filtered_SP = removed_sig_SP;
    
close all
pat_2_residue=abs(fft(sig_residue));
figure;plot(pat_1_residue);
title("Patient 1 residue")
figure;plot(pat_2_residue);
title("Patient 2 residue")
res_matrix_avg = mean([pat_1_residue,pat_2_residue])




    