% compute the initial condition for 96
% 2022.10.14 Kurama Okubo
clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)

%% 

cp = 6200; %[m/s]
cs = 3600; %[m/s]
rho = 2980; %[kg/m3]
h = 0.1;
T = 1e-7; % delta time

xdist_all = [0.05, 0.1, 0.15]; % horizontal distance from source to sensor (m)
zdist_all_bd = [h, 0, h]; % virtical dist from source to sensor (should be 0 or h) (m)
zdist_all_mt = [h/2, h/2, h/2]; % virtical dist from source to sensor (should be 0 or h) (m)
sensor_z = 0.07; % distance from the fault plane (m)

%%
fid = fopen('crossvalidation_distance_r.txt','w');
for sensorid = 1:3
    xdist = norm([xdist_all(sensorid), sensor_z], 2);
    zdist = zdist_all_mt(sensorid);
    fprintf(fid, "x:%12.8e z:%12.8e\n", xdist*1e-3, zdist*1e-3);
end
fclose(fid)
