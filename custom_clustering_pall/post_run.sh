#!/bin/bash

echo "Ensemble Clustering: Running Test Case 1"

echo "Downloading bootid files from pod"

scp -i ~/Documents/cloud_keys/cloud_related/POD/pod_rsa vinaya@192.41.74.245:/home/vinaya/Documents/MATLAB/Algorithm/Matlab/mytoolbox/clustering/methods_custom_classifiers/custom_clustering/tmp/bootId* /home/vinaya/Documents/matlab_custom_clust/tmp/ 

echo "Copying files from tmp/ --> unit_tests/tmp"

cp -R tmp/bootId* unit_tests/tmp/

exit 0
