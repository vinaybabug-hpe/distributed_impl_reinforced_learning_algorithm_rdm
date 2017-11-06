#!/bin/bash

echo "Ensemble Clustering: Running Test Case 1"

echo "Copying files from unit_tests/tmp --> tmp/"

cp -R unit_tests/tmp/bootId* tmp/

echo "Copying bootid files to pod"

scp -i ~/Documents/cloud_keys/cloud_related/POD/pod_rsa /home/vinaya/Documents/matlab_custom_clust/tmp/bootId* vinaya@192.41.74.245:/home/vinaya/Documents/MATLAB/Algorithm/Matlab/mytoolbox/clustering/methods_custom_classifiers/custom_clustering/tmp

echo "Copying *.sub files to pod"

scp -i ~/Documents/cloud_keys/cloud_related/POD/pod_rsa /home/vinaya/Documents/matlab_custom_clust/pod_*.sub vinaya@192.41.74.245:/home/vinaya/Documents/MATLAB/Algorithm/Matlab/mytoolbox/clustering/methods_custom_classifiers/custom_clustering/



exit 0
