addpath(genpath('/home/vinaya/Documents/MATLAB/Algorithm/Matlab'));


load('data/X_1000.mat');
load('data/bootIdxs_1000.mat');

load('data/PST.mat');


output_file=strcat(pwd,'/output_idx_1000.txt');


mex2c_cluster_fit_pall(X, PST, bootIdxs,output_file,48);

% C=mex2c_cluster_fit(X, PST, bootIdxs);
% B=struct();
% B.partitions=struct2dataset(C.partitions);
