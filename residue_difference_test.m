function [outputArg1,outputArg2] = residue_difference_test(varargin)
% Author: Alex Topping
%   Function shows the average and difference of two (todo: or more)
%   patient data files
%   Input a list of excel files in quotations


%%%% Bradley University %%% Yufeng Lu
%% revised for sEMG project 
%%% created: 9/7/2016
%%% revised:  5/24/2018 



%% prepare the data 

%%% initial:  when the subject came to the lab
%%% 10-min resting time: 10 mins after the initial capture
%% for initial or 10-min resting,  20 measurements are done. 10 for right side, 10 for left side
%% For example: R1_5:   right side, initial data capture,  the 5th data set


 %filename = 'init_test.mat';  %% need change the file name so that all results can be saved.
 
 %file = matfile(filename, 'Writable', true);

 paient_index = varargin;
 num_patients = numel(paient_index);
 
 
for patient_x = 1: num_patients
        %%%%%%%%% 
        f_sampling = 2048; %%sampling rate
        n_samples = f_sampling*5; %% pick 5 second long data
        xls_begin = 9; %% start at sheet B9
        xls_end = n_samples-1+xls_begin; %% select n_samples

        %string_xls_name = ['103-normal-M-6-6-13.xlsx']; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS

        string_xls_name = [paient_index(patient_x)]; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS


        %string_sheet = ['R1_5'];   % get the data from R1_5 sheet  

        string_sheet = ['R1_6'];   % get the data from R1_6 sheet  

        %string_sheet = ['R1_7'];   % get the data from R1_6 sheet  


        string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
        raw = -xlsread(string(string_xls_name),string_sheet,string_cells); %% 

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


         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  reconstructed_sig_EMD


        sig_seg = zeros(NEcho,NSamples);
        window_size = round(0.1 *f_sampling); %%% 0.1 second peak offset 

        sig_residue = reconstructed_sig_EMD; %%%%%%%%%%%%%%
        sig  = reconstructed_sig_EMD; %%%%%%%%%%%%%%

        peak_loc_seg = zeros(NEcho,1);


        STFT_M = 256;
        STFT_fm = f_sampling;


            [orig_T,orig_F,orig_SP] = STFT([t1; original_sig]', STFT_M, STFT_fm);

            file.orig_T = orig_T;
            file.orig_F = orig_F;
            file.orig_SP = orig_SP;


            [reconstructed_T,reconstructed_F,reconstructed_SP] = STFT([t1; reconstructed_sig_EMD]', STFT_M, STFT_fm);

            file.reconstructed_T = reconstructed_T;
            file.reconstructed_F = reconstructed_F;
            file.reconstructed_SP = reconstructed_SP;

            [residue_T,residue_F, residue_SP] = STFT([t1; sig_residue]', STFT_M, STFT_fm);

            file.filtered_T = residue_T;
            file.filtered_F = residue_F;
            file.filtered_SP = residue_SP;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            removed_sig = original_sig - reconstructed_sig_EMD; 

           [removed_sig_T,removed_sig_F, removed_sig_SP] = STFT([t1; removed_sig]', STFT_M, STFT_fm);

            file.filtered_T = removed_sig_T;
            file.filtered_F = removed_sig_F;
            file.filtered_SP = removed_sig_SP;
            
       figure;plot(sig_residue);
       title("placeholder title")
       
       text(x, y, sprintf('Text %f more text', variable))
            
            
end %ends the for loop
    


            
            
            
            
pat_1_residue = abs(fft(sig_residue));
pat_2_residue=abs(fft(sig_residue));


figure;plot(pat_1_residue);
title("Patient 1 residue")
figure;plot(pat_2_residue);
title("Patient 2 residue")
total_residue = pat_1_residue + pat_2_residue;
mean_residue = total_residue ./ 2; % 2 must be paramaterized for averages of more than 2 patients
figure;plot(mean_residue);
title("Average Residue of Patient 1 and Patient 2") % Verified on patient 116 and 225





    