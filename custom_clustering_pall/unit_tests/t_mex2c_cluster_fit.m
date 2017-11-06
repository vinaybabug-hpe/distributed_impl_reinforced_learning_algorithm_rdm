addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));

load('data/X_200.mat');
load('data/bootIdxs_200.mat');
% load('data/X_400.mat');
% load('data/bootIdxs_400.mat');
% load('data/X_600.mat');
% load('data/bootIdxs_600.mat');
% load('data/X_800.mat');
% load('data/bootIdxs_800.mat');
% load('data/X_1000.mat');
% load('data/bootIdxs_1000.mat');

load('data/PST.mat');

output_file=strcat(pwd,'/output_idx_200.txt');
%output_file=strcat(pwd,'/output_idx_400.txt');
%output_file=strcat(pwd,'/output_idx_600.txt');
%output_file=strcat(pwd,'/output_idx_800.txt');
%output_file=strcat(pwd,'/output_idx_1000.txt');


mex2c_cluster_fit_pspec(X, PST, bootIdxs,output_file,12);


