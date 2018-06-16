function [outputArg1,outputArg2] = avg_residue(pat1, pat2)
% Author: Alex Topping
%   Function shows the average and difference of two (todo: or more)
%   patient data files
%   Where x is the name patient 1 xcel data
%   Where y is the name of patient 2 xcel data


%%%% Bradley University %%% Yufeng Lu
%% revised for sEMG project 
%%% created: 9/7/2016
%%% revised:  5/24/2018 
temp1 = pat1
temp2 = pat2
%clear all; clc; close all;

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

string_xls_name = [temp1]; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS


%string_sheet = ['R1_5'];   % get the data from R1_5 sheet  

string_sheet = ['R1_6'];   % get the data from R1_6 sheet  

%string_sheet = ['R1_7'];   % get the data from R1_6 sheet  


string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
raw = -xlsread(string_xls_name,string_sheet,string_cells); %% 

figure;
plot(raw);grid on; title(['recorded data'  string_xls_name]);
xlabel('# sample ');
ylabel('Amplitude [mV]');




    
% Begin Data calculations for patient 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 
%%%%%%%%% 
f_sampling = 2048; %%sampling rate
n_samples = f_sampling*5; %% pick 5 second long data
xls_begin = 9; %% start at sheet B9
xls_end = n_samples-1+xls_begin; %% select n_samples

%string_xls_name = ['103-normal-M-6-6-13.xlsx']; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS

string_xls_name = [temp2]; %% file name: 103M-normal, 115F-normal, 203F-AS,208M-AS


%string_sheet = ['R1_5'];   % get the data from R1_5 sheet  

string_sheet = ['R1_6'];   % get the data from R1_6 sheet  

%string_sheet = ['R1_7'];   % get the data from R1_6 sheet  


string_cells = ['B' num2str(xls_begin) ':B' num2str(xls_end)]; % range of data cells
raw = -xlsread(string_xls_name,string_sheet,string_cells); %% 


figure;
plot(raw);grid on; title(['recorded data'  string_xls_name]);
xlabel('# sample ');
ylabel('Amplitude [mV]');



