% This is a Matlab scource code. Before executing it in Matlab, please 
% change the file extention from 'txt' to 'm', and put the 'polar_data.csv' file
% in the same folder.
clear
clc

aero_data = readtable('polar_data.csv', 'PreserveVariableNames' , true);

%% print max/min values

for i = 1:size(aero_data,2)
    col = aero_data.Properties.VariableNames(i);
    fprintf('The max value of  %s is: %#.3g\n', col{1}, max(aero_data{:,i}));
    fprintf('The min value of  %s is: %#.3g\n', col{1}, min(aero_data{:,i}));
end

%% compute and plot lift coefficients
avg = 5; % mean
amp = 3; % amplitude
T = 10; % period
fs = 10; % sampling frequency
time = 0:1/fs:T; 
angle_of_attack = avg + amp * sin(2*pi/T*time); 
lift_coefficient = interp1(aero_data{:,1}, aero_data{:,2}, ...
    angle_of_attack,'linear'); % interpolation based on suppied data

plot(time, lift_coefficient,'-o','MarkerSize', 4);
xlabel('Time (s)');
ylabel('Lift Coefficient');
grid on;