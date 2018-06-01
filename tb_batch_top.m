clear variables; close all; clc;
set(0,'DefaultFigureVisible', 'on');
set(0,'defaultfigurecolor',[1 1 1])

patientsToRun = [ 1, 14, 24, 28 ];

% 1 - female without AS (22)
% 14 - male without AS (18)
% 24 - female with AS (38)
% 28 - male with AS (32)


%for j = 41:1:47
for j = 2:1:2
    tb_semg_BigBatch( patientsToRun(j) );
    
    disp('completed: ');
    disp(j);
end

disp ('done!');