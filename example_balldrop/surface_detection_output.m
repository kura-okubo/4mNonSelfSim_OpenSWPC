% Script to debug the index of surfaces in OpenSWPC
% to assign the second-order computation
% 2022/06/27 Kurama Okubo

% 2025/03/12 update for master plot

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

%% simulateion parameters for visualization
dx = 2.0e-3; % [m] grid size
dy = 2.0e-3; % [m] grid size
dz = 2.0e-3; % [m] grid size

%% read debug output of surfaces
rootdir = "./out";

kfs_biax_upper = importdata(rootdir+"/kfs_biax_upper.dat");
kfs_biax_lower = importdata(rootdir+"/kfs_biax_lower.dat");
kfs_biax_upper_top = importdata(rootdir+"/kfs_biax_upper_top.dat");
kfs_biax_upper_bot = importdata(rootdir+"/kfs_biax_upper_bot.dat");
kfs_biax_lower_top = importdata(rootdir+"/kfs_biax_lower_top.dat");
kfs_biax_lower_bot = importdata(rootdir+"/kfs_biax_lower_bot.dat");

%% Plot k indices of surface
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 800 500];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on;

plot3(kfs_biax_upper.data(:, 4)*dx, kfs_biax_upper.data(:, 5)*dy, kfs_biax_upper.data(:, 6)*dz, '.', "Color", "b", "DisplayName","kfs\_biax\_upper", "MarkerSize",12);
plot3(kfs_biax_lower.data(:, 4)*dx, kfs_biax_lower.data(:, 5)*dy, kfs_biax_lower.data(:, 6)*dz, '.', "Color", "r", "DisplayName","kfs\_biax\_lower", "MarkerSize",12);
plot3(kfs_biax_upper_top.data(:, 4)*dx, kfs_biax_upper_top.data(:, 5)*dy, kfs_biax_upper_top.data(:, 6)*dz, 's', "Color", [0, 0.7, 0], "DisplayName","kfs\_biax\_upper\_top", "MarkerSize",12);
plot3(kfs_biax_upper_bot.data(:, 4)*dx, kfs_biax_upper_bot.data(:, 5)*dy, kfs_biax_upper_bot.data(:, 6)*dz, '.', "Color", 'c', "DisplayName","kfs\_biax\_upper\_bot", "MarkerSize",12);
plot3(kfs_biax_lower_top.data(:, 4)*dx, kfs_biax_lower_top.data(:, 5)*dy, kfs_biax_lower_top.data(:, 6)*dz, '.', "Color", [0.8, 0.5, 0.4], "DisplayName","kfs\_biax\_lower\_top", "MarkerSize",12);
plot3(kfs_biax_lower_bot.data(:, 4)*dx, kfs_biax_lower_bot.data(:, 5)*dy, kfs_biax_lower_bot.data(:, 6)*dz, '.', "Color", 'y', "DisplayName","kfs\_biax\_lower\_bot", "MarkerSize",12);

% 
% xlim([140, 160]);
% xlim([-10, 10]);
% ylim([-100, 100]);
% zlim([-20, 120]);
xlabel("x [m]");
ylabel("y [m]");
zlabel("z [m]");

legend("Location","best");
view([121 20]);
set(gca, 'Zdir', 'reverse')

set(gcf, 'Color', 'w');
figname = "./kfs_biax.png";
exportgraphics(fig, figname,  'Resolution', 100);
