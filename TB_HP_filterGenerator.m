function [ fil_Num,fil_Dem ] = TB_HP_filterGenerator( fc, tb, fs, Rs, Rp, lengthH, choice)

%% design a HP filter to apply in the time domain and return the impulse response
% also, plot the filter response, normalized against the nyquist freq.
%% VARIABLE SUMMARY
% fcs = cutoff freq [Hz]
% tb = transition band width [Hz]
% fs = sampling freq. [Hz]
% Rs = Stopband Attenuation Ripple Max [dB]
% Rp = Passband Attenuation Ripple Max [dB]
% lengthH = length of impulse response [#samples]
% Choice picks the filter type [integer]


plotsOn = 1;
%% Begin Function
% cutoff frequencies
fcs = fc-tb;
fcp = fc;  

% normalize against Nyquist freq
Wp = fcp*2/fs;
Ws = fcs*2/fs;

switch choice
   case 1
%% Butterworth digital filter design
    type = 'Butterworth,';
    [filt_Order,Wn_butt] = buttord(Wp,Ws,Rp,Rs);
    [fil_Num,fil_Dem] = butter(filt_Order,Wn_butt,'high');   
   case 2
%% Chebyshev type I filter
    type = 'Chebyshev Type I,';
    [filt_Order,Wn_cheb1] = cheb1ord(Wp,Ws,Rp,Rs);
    [fil_Num,fil_Dem] = cheby1(filt_Order,Rp,Wp,'high'); 
   
   case 3
       type = 'Chebyshev Type II,';
    [filt_Order,Wn_cheb2] = cheb2ord(Wp,Ws,Rp,Rs);
    [fil_Num,fil_Dem] = cheby2(filt_Order,Rs,Ws,'high');

   case 4
%% Elliptic filter
    type = 'Elliptic,';
    [filt_Order,Wn_ellip] = ellipord(Wp,Ws,Rp,Rs);
    [fil_Num,fil_Dem] = ellip(filt_Order,Rp,Rs,Wp,'high');
   otherwise
      
end


% need freqz for the plot
[H_filt,W_filt]= freqz(fil_Num,fil_Dem,fs);
% impz gets the impulse response for use in STFT
[impRes,t1] = impz(fil_Num,fil_Dem, lengthH);
 
    if(plotsOn  == 1)
        figure
        subplot(2,1,1)
        plot(W_filt*fs/(2*pi),abs(H_filt),'linewidth',2,'color',[0,0.7,0])
        grid on, box on
        title(['Digital IIR Filter design: ' type ' Order: ' num2str(filt_Order) ' (Amplitude Response)']);
        xlabel('Normalized Frequency(\times\pi rad/sample)');
        ylabel('Amplitude');
        ylim([0 1.2]);
        set(gca,'XTick',[0 Wp 1]*fs/2)
        set(gca,'ytick',[0 1 1+10^(Rp/(-20))])
        %text(fcp,1-(10^(Rp/(-20))+0.01),['Passband ripple = ',num2str(10^(Rp/(-20)))]);
        plot_settings( [0 0 0] )
        
        subplot(2,1,2)
        plot(t1/fs, impRes);
        grid on, box on
        title(['Digital IIR HP Filter: ' type ' Order: ' num2str(filt_Order) ' (Imp. Response)']);
        xlabel('Time [sec]');
        ylabel('Amplitude');
        plot_settings( [0 0 0] )
        
    end
end

