function [outputArg1,outputArg2] = avgmax()
%	Author: Alex Topping
%   Plots multiple patients and finds averages as well as max points


%% Plotting individual patients

% Patient 1
load("patient 103M-AS sheet R2_5.mat", "sig_residue") 
close
res1 = abs(fft(sig_residue, 2048));
res1norm = res1 ./ max(res1); 

% Patient 2
load("patient 103M-AS sheet R2_4.mat", "sig_residue") 
close
res2 = abs(fft(sig_residue, 2048));
res2norm = res2 ./ max(res2); 

% Patient 3
load("patient 103M-AS sheet R2_3.mat", "sig_residue") 
close
res3 = abs(fft(sig_residue, 2048));
res3norm = res3 ./ max(res3); 

% Patient 4
load("patient 103M-AS sheet R2_2.mat", "sig_residue") 
close
res4 = abs(fft(sig_residue, 2048));
res4norm = res4 ./ max(res4); 

% Patient 5
load("patient 103M-AS sheet R1_5.mat", "sig_residue") 
close
res5 = abs(fft(sig_residue, 2048));
res5norm = res5 ./ max(res5); 

% Patient 6
load("patient 103M-AS sheet R1_4.mat", "sig_residue") 
close
res6 = abs(fft(sig_residue, 2048));
res6norm = res6 ./ max(res6); 

% Patient 7
load("patient 103M-AS sheet R1_3.mat", "sig_residue") 
close
res7 = abs(fft(sig_residue, 2048));
res7norm = res7 ./ max(res7); 

% Patient 8
load("patient 103M-AS sheet R1_2.mat", "sig_residue") 
close
res8 = abs(fft(sig_residue, 2048));
res8norm = res8 ./ max(res8); 

% Patient 9
load("patient 103M-AS sheet R1_1.mat", "sig_residue") 
close
res9 = abs(fft(sig_residue, 2048));
res9norm = res9 ./ max(res9); 

% Patient 10
load("patient 103M-AS sheet L2_10.mat", "sig_residue") 
close
res10 = abs(fft(sig_residue, 2048));
res10norm = res10 ./ max(res10); 

% Patient 11
load("patient 103M-AS sheet L2_9.mat", "sig_residue") 
close
res11 = abs(fft(sig_residue, 2048));
res11norm = res11 ./ max(res11); 

% Patient 12
load("patient 103M-AS sheet L2_8.mat", "sig_residue") 
close
res12 = abs(fft(sig_residue, 2048));
res12norm = res12 ./ max(res12); 

% Patient 13
load("patient 103M-AS sheet L2_7.mat", "sig_residue") 
close
res13 = abs(fft(sig_residue, 2048));
res13norm = res13 ./ max(res13); 

% Patient 14
load("patient 103M-AS sheet L2_6.mat", "sig_residue") 
close
res14 = abs(fft(sig_residue, 2048));
res14norm = res14 ./ max(res14); 

% Patient 15
load("patient 103M-AS sheet L2_5.mat", "sig_residue") 
close
res15 = abs(fft(sig_residue, 2048));
res15norm = res15 ./ max(res15); 

% Patient 16
load("patient 103M-AS sheet L2_4.mat", "sig_residue") 
close
res16 = abs(fft(sig_residue, 2048));
res16norm = res16 ./ max(res16); 

% Patient 17
load("patient 103M-AS sheet L2_3.mat", "sig_residue") 
close
res17 = abs(fft(sig_residue, 2048));
res17norm = res17 ./ max(res17); 

% Patient 18
load("patient 103M-AS sheet L2_2.mat", "sig_residue") 
close
res18 = abs(fft(sig_residue, 2048));
res18norm = res18 ./ max(res18); 

% Patient 19
load("patient 103M-AS sheet L2_1.mat", "sig_residue") 
close
res19 = abs(fft(sig_residue, 2048));
res19norm = res19 ./ max(res19); 

% Patient 20
load("patient 103M-AS sheet L1_10.mat", "sig_residue") 
close
res20 = abs(fft(sig_residue, 2048));
res20norm = res20 ./ max(res20); 

% Patient 21
load("patient 103M-AS sheet L1_9.mat", "sig_residue") 
close
res21 = abs(fft(sig_residue, 2048));
res21norm = res21 ./ max(res21); 

% Patient 22
load("patient 103M-AS sheet L1_8.mat", "sig_residue") 
close
res22 = abs(fft(sig_residue, 2048));
res22norm = res22 ./ max(res22); 

% Patient 23
load("patient 103M-AS sheet L1_7.mat", "sig_residue") 
close
res23 = abs(fft(sig_residue, 2048));
res23norm = res23 ./ max(res23); 

% Patient 24
load("patient 103M-AS sheet L1_6.mat", "sig_residue") 
close
res24 = abs(fft(sig_residue, 2048));
res24norm = res24 ./ max(res24); 

% Patient 25
load("patient 103M-AS sheet L1_5.mat", "sig_residue") 
close
res25 = abs(fft(sig_residue, 2048));
res25norm = res25 ./ max(res25); 

% Patient 26
load("patient 103M-AS sheet L1_4.mat", "sig_residue") 
close
res26 = abs(fft(sig_residue, 2048));
res26norm = res26 ./ max(res26); 

% Patient 27
load("patient 103M-AS sheet L1_3.mat", "sig_residue") 
close
res27 = abs(fft(sig_residue, 2048));
res27norm = res27 ./ max(res27); 

% Patient 28
load("patient 103M-AS sheet L1_2.mat", "sig_residue") 
close
res28 = abs(fft(sig_residue, 2048));
res28norm = res28 ./ max(res28); 

% Patient 29
load("patient 103M-AS sheet L1_1.mat", "sig_residue") 
close
res29 = abs(fft(sig_residue, 2048));
res29norm = res29 ./ max(res29); 

% Patient 30
load("patient 103M-AS sheet R2_1.mat", "sig_residue") 
close
res30 = abs(fft(sig_residue, 2048));
res30norm = res30 ./ max(res30); 


 figure;plot(res1norm)
 figure;plot(res2norm) 
% figure;plot(res3norm) 
% figure;plot(res4norm) 
% figure;plot(res5norm) 
% figure;plot(res6norm) 
% figure;plot(res7norm) 
% figure;plot(res8norm) 
% figure;plot(res9norm) 
% figure;plot(res10norm) 
% figure;plot(res11norm) 
% figure;plot(res12norm) 
% figure;plot(res13norm) 
% figure;plot(res14norm) 
% figure;plot(res15norm)
% figure;plot(res16norm)
% figure;plot(res17norm) 
% figure;plot(res18norm) 
% figure;plot(res19norm) 
% figure;plot(res20norm) 
% figure;plot(res21norm) 
% figure;plot(res22norm) 
% figure;plot(res23norm) 
% figure;plot(res24norm) 
% figure;plot(res25norm) 
% figure;plot(res26norm) 
% figure;plot(res27norm) 
% figure;plot(res28norm) 
% figure;plot(res29norm) 
% figure;plot(res30norm) 


%% Calculate an average
total_res = res1norm + res2norm + res3norm + res4norm + res5norm + res6norm + res7norm + res8norm + res9norm + res10norm + res11norm + res12norm + res13norm + res14norm + res15norm + res16norm + res17norm + res18norm + res19norm + res20norm + res21norm + res22norm + res23norm + res24norm + res25norm + res26norm + res27norm + res28norm + res29norm + res30norm; % Add up 3 signals
mean_res = total_res ./ 30; % Divide by 30 since there are 30 signals
mean_res_norm = mean_res ./ max(mean_res);


%% Plot max points
[pks, locs] = findpeaks(mean_res);
[sort_peak, sort_index] = sort(pks, 'descend');

max5freq = locs(sort_index(1:5)) .* 2 % Multiply by 2 because of FFT index


figure;plot(mean_res_norm)
title("Average")




end

