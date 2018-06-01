function [T,F,SP] = STFT(D, M, fm)


%% Revised STFT
%% Y. Lu  7/20/2014
% y: input signal
% D: the data matrix: D(:,1)=t; D(:,2)=y;
% M: the lenght of the FT window
% [0 fm] The Frequency range for the signal
% SP: the spectogram (short time fourier transform) of D 

t=D(:,1);y=D(:,2);
N=length(y);

W=hanning(M);
ye=zeros(N+M,1);
ye(M/2+1:M/2+N)=y;
%yea=hilbert(ye);
dt=t(2)-t(1);% time resolution
fs=1/dt; % the sampling frequency

df=fs/M; %frequency resolution
L=round(fm/df); % the index of the max freq.
f=[0:L-1]*df;

i=1;
while i<=N
   
   s=ye(i:i+M-1).*W;
   sa=hilbert(s);
   
   S=fft(sa);
   B(:,i)=S(1:L);
   i=i+1;% increase by one (next sample)
end

SP = (1/N)*abs(B.*B);
SP = SP/max(max(SP)); % normalize it
T  = t;
F  = f;

% mesh(T,F,SP)
% contour(T,F,SP)
