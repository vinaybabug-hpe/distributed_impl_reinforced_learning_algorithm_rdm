addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));

load('data/X_200.mat');
load('data/bootIdxs_200.mat');


load('data/PST.mat');

output_file=strcat(pwd,'/output_idx_200.txt');



mex2c_cluster_fit_pspec(X, PST, bootIdxs,output_file,8);

% C=mex2c_cluster_fit(X, PST, bootIdxs);
% B=struct();
% B.partitions=struct2dataset(C.partitions);
