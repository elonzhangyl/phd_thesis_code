addpath 'C:\loong\MATLAB WORKPLACE\6.distribution_linear_no_water_real_new_noise\area3'
[X,Y]=meshgrid(1:30,1:5);
options.barwidth = 0.3;
% options.TScaling=[0.9 0.9;0.9 0.9;0.9 0.9];
area3(X,1:5,squeeze(avg_AUC(1,:,:))',options)

[X,Y]=meshgrid(1:5,1:30);
area3(X,1:30,squeeze(avg_AUC(1,:,:)))

size(X)
size(squeeze(avg_AUC(1,:,:)))