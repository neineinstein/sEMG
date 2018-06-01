function [ ] = tb_acquirePatientData( age, PatientID, AS_info, sex, string_xls_name)

%for CFSweep=20:20:60 % iterate through this function with different cutoff frequencies
    
%close all;   % computer memory saver   

color1 = [128 0 0]/255;
color2 = [0   0 128]/255;
color3 = [190 190 190]/255;

fs = 2048;          % sampling frequency of instrument
STFT_M = 128;
STFT_FM = 1024;
% settings for the time domain filter applied to raw signal
low_CutoffFrequency = 90;           % cutoff freq [Hz]
transitionBandFreqWidth = 10;            % width of the transition band [Hz] (smaller values here make a longer filter)
stopbandRippleAttenuation = 80;            % stop band attentuation [dB down]
passbandAttenuationRipple = 0.01;          %pass band attenuation ripple (dB flex about 1)
filterChoice = 4;         % 1 - butter, 2 - cheby1, 3 - cheby2, 4 - elliptic
tolerance = 0.2;          % [%] of max value used in peak finding
Imp_Length = 250;    % biggest fudge factor in project...length in samples of the filter

filepath = './data';
printFileName = ['(P',PatientID,')(A-',age,')(S-',sex,')(AS-',AS_info,')(CFHz-',num2str(low_CutoffFrequency),')'];
% fileID = fopen([filepath,'/Pat_',num2str(PatientID),'/FreqReadOut',printFileName,'.csv'],'w');
fileID = fopen([filepath,'FreqReadOut',printFileName,'.csv'],'w');

   %% generate a filter        
[ fil_Num,fil_Dem ] = TB_HP_filterGenerator( low_CutoffFrequency, transitionBandFreqWidth, fs, stopbandRippleAttenuation, passbandAttenuationRipple, Imp_Length, filterChoice);

[ sheetRange ] = getSheet();

%% lets try to make a cell to hold all of the FFTs that will be generated.  
%  The cell needs to 
% 1.) hold the sheet,
% 2.) the TD filtered segments, and the
% 3.) segments
% so the depth of the cell is the length of the sheetRange index
[~,cellDepth] = size(sheetRange);
maxNumberEventsInSignal = 10;

ecg_width = 500;
qz_width = 500;   
f_qz = fs*(0:floor((qz_width/2)))/qz_width;
f_ecg = fs*(0:floor((ecg_width/2)))/ecg_width;
filtered_QZ_FFTs = zeros(maxNumberEventsInSignal, length(f_qz), cellDepth);
raw_QZ_FFTs = zeros(maxNumberEventsInSignal, length(f_qz), cellDepth);
filtered_ECG_FFTs = zeros(maxNumberEventsInSignal, length(f_ecg), cellDepth);
raw_ECG_FFTs = zeros(maxNumberEventsInSignal, length(f_ecg), cellDepth);


for loopIndex = 1:1:length(sheetRange)
    
    n_samples = fs*5;   % length of the segment used in the analysis
    xls_begin = 9;      % start at B9 on excel sheet, double check input from here if it looks weird
    xls_end = n_samples-1+xls_begin; %% select n_samples
    dataTrial = char(sheetRange(loopIndex));
    string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
    raw =xlsread(string_xls_name,dataTrial,string_cells); % generate input data

    %[ s ] = tb_HP_Filter_Struct(low_CutoffFrequency,transitionBandFreqWidth,stopbandRippleAttenuation,passbandAttenuationRipple, filterChoice, tolerance, Imp_Length );
    
    filepathIter = [filepath,strcat('/Pat_',num2str(PatientID),'/',dataTrial,'-CFHz_',num2str(low_CutoffFrequency),'-STFTM_',num2str(STFT_M))];
    mkdir(filepathIter);
    
    samples = (1:1:length(raw));
    [ecgIndex, qzIndex] = tb_investigateZone( raw, samples, fs, ecg_width, qz_width);

    ecgZ_TD = zeros(length(ecgIndex),ecg_width+1);
    [ ecgZ, ecgZX ] = zoneBuilder(ecgIndex,ecg_width, raw, samples);
   
    qzZ_TD = zeros(length(qzIndex),qz_width+1);    
    [ qZ, qZX ] = zoneBuilder(qzIndex,qz_width, raw, samples);

    [m, n] = size(qZ);
    
    % this loop is filtering out bad QZ acquisitions by zeroing them out.
    % later functions look for all zero data rows and discard them.
    flagData = 0;
    for j=1:1:m
        for k=1:1:n
            if(qZ(j,k) >= 20)
                flagData = 1;
            end
        end
        %zero out this data line
         for k=1:1:n
             if(flagData == 1)
                qZ(j,k) = 0;
            end
         end
        flagData = 0;
    end

   %% generate time domain filtered segments
    for j = 1:1:length(ecgIndex)
             ecgZ_TD(j,:) = filter(fil_Num,fil_Dem, ecgZ(j,:));
    end

    for j = 1:1:length(qzIndex)
          qzZ_TD(j,:) = filter(fil_Num,fil_Dem, qZ(j,:));
    end 
    
%        segName = filepathIter;
%        mkdir(segName);
%        disp('View Segmentation Results');
        %% plot the segmentations
            figure
            plot(samples/fs, raw, 'Color',color3,'linewidth',0.5);  hold on;
            plot(qZX(1,:)/fs,qZ(1,:),'Color',color2,'linewidth',2); 
            plot(ecgZX(1,:)/fs,ecgZ(1,:),'Color',color1,'linewidth',2);
            legend('Unused','Quiet Zones','Loud Zones');
            for j = 2:1:length(qzIndex)
                plot(qZX(j,:)/fs,qZ(j,:),'Color',color2,'linewidth',2); hold on;
            end

            for j = 2:1:length(ecgIndex)
                plot(ecgZX(j,:)/fs,ecgZ(j,:),'Color',color1,'linewidth',2); hold on;
            end
            hold off;
%             title(['(Patient',PatientID,')(Side-',dataTrial(1),')(Sex-',sex,')(AS-',AS_info,')']);
            xlabel('Time [sec]');
            ylabel('Amplitude [mV]');
            %ax = gcf;
            plot_settings( [0, 0, 0] );
%            tb_plotSaver( gcf, [printFileName,'SegResults'], segName )  
       % close all;
    
       
        disp(['(Patient',PatientID,')(Side-',dataTrial(1),')(Sex-',sex,')(AS-',AS_info,')']);
       disp('View STFT Results - Individual');
% turned these off for the latest part of the code.  They
% add alot of processing time
        segName = filepathIter;
        mkdir(segName);
        % yLo/yHi are frequencies used in plotting windows, they have no mathematical effect
        yLo = 0;
        yHi = 400;  
        yHi = fs/2;
            for j=1:1:length(ecgIndex)
                 plotSTFT(STFT_M, STFT_FM, ecgZ(j,:), ecgZ_TD(j,:), ecgZX(j,:), fs, low_CutoffFrequency ,j, 'LZ', yLo, yHi, printFileName, segName);
            end

%             close all;

        yLo = 0;
        yHi = 400;
       yHi = fs/2;
            for j=1:1:length(qzIndex)
                 plotSTFT(STFT_M, STFT_FM, qZ(j,:), qzZ_TD(j,:), qZX(j,:), fs,  low_CutoffFrequency ,j, 'QZ', yLo, yHi, printFileName, segName);
            end

%        close all;

       %filepathIter = [filepath,strcat('/Pat_',num2str(PatientID),'/',dataTrial,'/CFHz_',num2str(low_CutoffFrequency))];

%        disp('View FFT Results - Individual'); 
       segName = filepathIter;
       mkdir(segName);
%        fprintf(fileID,[dataTrial,'\n']);
       [ ecgFFT_tdFilt, ecgFFT, ~ ] = fftManager( ecgZ, ecgZ_TD, fs, tolerance, 'LZ',printFileName, segName, dataTrial, fileID );
       [ qzFFT_tdFilt, qzFFT, ~ ] = fftManager( qZ, qzZ_TD, fs, tolerance, 'QZ',printFileName, segName, dataTrial, fileID );
       
       [szECGy, szECGx] = size(ecgFFT);
       [szQZy, szQZx] = size(qzFFT);
       
       for k=1:1:szECGy
          filtered_ECG_FFTs(k, :, loopIndex) = ecgFFT_tdFilt(k,:);
          raw_ECG_FFTs(k, :, loopIndex) = ecgFFT(k,:);
       end
       
       for k=1:1:szQZy
          filtered_QZ_FFTs(k, :, loopIndex) =  qzFFT_tdFilt(k,:);
          raw_QZ_FFTs(k, :, loopIndex) =  qzFFT(k,:);
       end
%        close all;
       
%end

       segName = filepathIter;
       mkdir(segName);
       % this function handles the frequency domain plots from the FFT
       plotFFTs( f_qz, f_ecg, raw_ECG_FFTs, filtered_ECG_FFTs, raw_QZ_FFTs, filtered_QZ_FFTs , printFileName, segName);


      

 fclose(fileID);
disp('done');
 
end