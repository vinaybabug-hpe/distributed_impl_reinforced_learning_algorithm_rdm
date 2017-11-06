/*
 ============================================================================
 Name        : kmeans_nwu.c
 Author      : Vinay Gavirangaswamy
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int _debug;
//#include "nwu/kmeans.h"
#include "univ_of_tokyo/cluster.h"


int main(void) {

	int i, j, isBinaryFile, is_output_timing, verbose;

	int numClusters, numCoords, numObjs;
	float **objects;  /* [numObjs][numCoords] data objects */
	int *idx;
	float **clusters; /* [numClusters][numCoords] cluster center */
	float threshold;
	char *filename, *center_filename;
	int    *membership;    /* [numObjs] */

	/* some default values */
	_debug = 0;
	verbose = 1;
	threshold = 0.001;
	numClusters = 5;
	isBinaryFile = 0;
	is_output_timing = 0;
	filename = "color100.txt";
	center_filename = NULL;



	objects = file_read(isBinaryFile, filename, &numObjs, &numCoords);


//	idx = (int*) malloc(numObjs * sizeof(int));

	/* allocate a 2D space for clusters[] (coordinates of cluster centers)
	 this array should be the same across all processes                  */
	clusters = (float**) malloc(numClusters * sizeof(float*));
	assert(clusters != NULL);
	clusters[0] = (float*) malloc(numClusters * numCoords * sizeof(float));

	assert(clusters[0] != NULL);
	for (i = 1; i < numClusters; i++)
		clusters[i] = clusters[i - 1] + numCoords;

	printf("selecting the first %d elements as initial centers\n", numClusters);
	/* copy the first numClusters elements in feature[] */
	for (i = 0; i < numClusters; i++)
		for (j = 0; j < numCoords; j++)
			clusters[i][j] = objects[i][j];
//
//	/* check initial cluster centers for repeatition */
//	if (check_repeated_clusters(numClusters, numCoords, clusters) == 0) {
//		printf(	"Error: some initial clusters are repeated. Please select distinct initial centers\n");
//		return 1;
//	}
//
//    /* start the timer for the core computation -----------------------------*/
//    /* membership: the cluster id for each data object */
//    membership = (int*) malloc(numObjs * sizeof(int));
//    assert(membership != NULL);
//
//    seq_kmeans(objects, numCoords, numObjs, numClusters, threshold, membership,
//               clusters);
//
//    free(objects[0]);
//    free(objects);
//
///*
//    printf("\nCluster Centroids\n");
//    for(i=0; i< numClusters; i++){
//    	for(j=0; j<numCoords; j++){
//    		printf(" %f ", clusters[i][j]);
//    	}
//
//    	printf("\n");
//    }
//*/
//
//
//    center_centroids_idx2(numClusters, numCoords, clusters, numObjs, numCoords, objects, idx);
//
////    printf("\nCluster IDX\n");
////
////    for (i = 0; i < numObjs; i++) {
////
////			printf(" %d ", idx[i]);
////
////
////		printf("\n");
////	}


	// University of Tokyo, University of Tokyo, Institute of Medical Science,
	// Human Genome Center, Laboratory of DNA Information Analysis library

	int** mask;
	double* weight = malloc(numCoords*sizeof(double));
	const int transpose = 0;
	const char dist = 'e';
	const char method = 'a';
	int npass = 1;
	int ifound = 0;
	double error;

	int* clusterid = malloc(numObjs*sizeof(int));

	for (i = 0; i < numCoords; i++) weight[i] = 1.0;



	/*kcluster(numClusters, numObjs, numCoords, objects, mask,weight,transpose,npass,method,dist,
		    clusterid, &error, &ifound);*/

	/*kmeans_adapter(numClusters, numObjs, numCoords, objects,
			     npass, method, dist, clusterid);*/

	/*agglom_adapter(numClusters, numObjs, numCoords, objects, method, dist, clusterid);*/

	kmedoids_adapter(numClusters,npass, numObjs, numCoords, objects, method, dist, clusterid);

	printf ("Cluster assignments:\n");
	for (i = 0; i < numObjs; i++)
	    printf ("Gene %d: cluster %d\n", i, clusterid[i]);

	return EXIT_SUCCESS;
}
