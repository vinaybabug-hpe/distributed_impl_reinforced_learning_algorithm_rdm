addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));


load('data/X_800.mat');
load('data/bootIdxs_800.mat');


load('data/PST.mat');


output_file=strcat(pwd,'/output_idx_800.txt');



mex2c_cluster_fit_pspec(X, PST, bootIdxs,output_file,8);


