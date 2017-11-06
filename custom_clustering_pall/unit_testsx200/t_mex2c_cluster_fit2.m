addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));


load('data/X2.mat');
load('data/PST2.mat');
load('data/bootIdxs2.mat');
output_file=strcat(pwd,'/output_idx.txt');
mex2c_cluster_fit_preEnsemble(X, PST, bootIdxs,output_file);

% C=mex2c_cluster_fit(X, PST, bootIdxs);
% B=struct();
% B.partitions=struct2dataset(C.partitions);
