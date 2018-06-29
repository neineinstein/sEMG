function [outputArg1,outputArg2] = test2()
%	Author: Alex Topping
%   Plots multiple patients and finds averages as well as max points


%% Plotting individual patients

% Patient 1
load("patient 103M-AS sheet L1_9.mat", "sig_residue") 
close
res1 = abs(fft(sig_residue, 2048));
res1norm = res1 ./ max(res1); 

% Patient 2
load("patient 103M-AS sheet L2_1.mat", "sig_residue") 
close
res2 = abs(fft(sig_residue, 2048));
res2norm = res2 ./ max(res2); 

% Patient 3
load("patient 103M-AS sheet R2_2.mat", "sig_residue") 
close
res3 = abs(fft(sig_residue, 2048));
res3norm = res3 ./ max(res3); 

% Patient 4
load("patient 103M-AS sheet L2_3.mat", "sig_residue") 
close
res4 = abs(fft(sig_residue, 2048));
res4norm = res4 ./ max(res4); 

% Patient 5
load("patient 103M-AS sheet L2_4.mat", "sig_residue") 
close
res5 = abs(fft(sig_residue, 2048));
res5norm = res5 ./ max(res5); 



 figure;plot(res1norm)
 figure;plot(res2norm) 
 figure;plot(res3norm) 
 figure;plot(res4norm) 
 figure;plot(res5norm) 


%% Calculate an average
total_res = res1norm + res2norm + res3norm + res4norm + res5norm; % Add up 5 signals
mean_res = total_res ./ 5; % Divide by 5 since there are 5 signals
mean_res_norm = mean_res ./ max(mean_res);


%% Plot max points
[pks, locs] = findpeaks(mean_res);
[sort_peak, sort_index] = sort(pks, 'descend');

max5freq = locs(sort_index(1:15)) .* 4 % Multiply by 2 because of FFT index


figure;plot(mean_res_norm)
title("Average")




end

