addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));

load('data/X_100.mat');
load('data/bootIdxs_100.mat');


load('data/PST.mat');

output_file=strcat(pwd,'/output_idx_100.txt');



mex2c_cluster_fit_pspec(X, PST, bootIdxs,output_file,1);


