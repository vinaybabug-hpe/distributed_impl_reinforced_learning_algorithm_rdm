/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_kmeans.c  (sequential version)                        */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*                                                                           */
/*   Copyright (C) 2005, Northwestern University                             */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nwu/kmeans.h"


/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
float euclid_dist_2(int    numdims,  /* no. dimensions */
                    float *coord1,   /* [numdims] */
                    float *coord2)   /* [numdims] */
{
    int i;
    float ans=0.0;

    for (i=0; i<numdims; i++)
        ans += (coord1[i]-coord2[i]) * (coord1[i]-coord2[i]);

    return(ans);
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         float  *object,      /* [numCoords] */
                         float **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    float dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}

/*----< seq_kmeans() >-------------------------------------------------------*/

/* Adapter return an array of cluster centers of size [numClusters][numCoords]       */
int seq_kmeans_adapter(double *objects,      /* in: [numObjs][numCoords] */
               int     numCoords,    /* no. features */
               int     numObjs,      /* no. objects */
               int     numClusters,  /* no. clusters */
               float   threshold,    /* % objects change membership */
               int *idx)     /* out: [numClusters][numCoords] */
{
	int    *membership;    /* [numObjs] */
	float **clusters, **x;
	int i,j;

    /* start the timer for the core computation -----------------------------*/
    /* membership: the cluster id for each data object */
    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);


	/* allocate a 2D space for clusters[] (coordinates of cluster centers)
	 this array should be the same across all processes                  */
	clusters = (float**) malloc(numClusters * sizeof(float*));
	clusters[0] = (float*) malloc(numClusters * numCoords * sizeof(float));
	assert(clusters[0] != NULL);
	for (i = 1; i < numClusters; i++)
		clusters[i] = clusters[i - 1] + numCoords;

	/* allocate a 2D space for data as x[] (coordinates of cluster centers)
		 this array should be the same across all processes                  */
	x = (float**) malloc(numObjs * sizeof(float*));
	x[0] = (float*) malloc(numObjs * numCoords * sizeof(float));
	for (i=1; i<numObjs; i++)
	            x[i] = x[i-1] + (numCoords);

	for(i=0; i<numObjs; i++){
		for(j=0;j<numCoords; j++){

			x[i][j] = objects[j+i*(int)numCoords];

		}
	}

//	printf("selecting the first %d elements as initial centers\n", numClusters);
		/* copy the first numClusters elements in feature[] */
		for (i = 0; i < numClusters; i++)
			for (j = 0; j < numCoords; j++)
				clusters[i][j] = x[i][j];

	seq_kmeans(x, numCoords, numObjs, numClusters, threshold, membership,
			clusters);

//    printf("\nCluster Centroids\n");
//    for(i=0; i< numClusters; i++){
//    	for(j=0; j<numCoords; j++){
//    		printf(" %f ", clusters[i][j]);
//    	}
//
//    	printf("\n");
//    }

	center_centroids_idx2(numClusters, numCoords, clusters, numObjs, numCoords,
			x, idx);

//	printf("\nCluster IDX\n");
//
//	for (i = 0; i < numObjs; i++) {
//
//		printf(" %d ", idx[i]);
//
//		printf("\n");
//	}

}

/* return an array of cluster centers of size [numClusters][numCoords]       */
int seq_kmeans(float **objects,      /* in: [numObjs][numCoords] */
               int     numCoords,    /* no. features */
               int     numObjs,      /* no. objects */
               int     numClusters,  /* no. clusters */
               float   threshold,    /* % objects change membership */
               int    *membership,   /* out: [numObjs] */
               float **clusters)     /* out: [numClusters][numCoords] */

{
    int      i, j, index, loop=0;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    float    delta;          /* % of objects change their clusters */
    float  **newClusters;    /* [numClusters][numCoords] */



    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    /* need to initialize newClusterSize and newClusters[0] to all 0 */
    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);

    newClusters    = (float**) malloc(numClusters * sizeof(float*));
    assert(newClusters != NULL);
    newClusters[0] = (float*)  calloc(numClusters * numCoords, sizeof(float));
    assert(newClusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        newClusters[i] = newClusters[i-1] + numCoords;
//



    do {
        delta = 0.0;

        for (i=0; i<numObjs; i++) {
            /* find the array index of nestest cluster center */
            index = find_nearest_cluster(numClusters, numCoords, objects[i],
                                         clusters);

            /* if membership changes, increase delta by 1 */
            if (membership[i] != index) delta += 1.0;

            /* assign the membership to object i */
            membership[i] = index;

            /* update new cluster center : sum of objects located within */
            newClusterSize[index]++;
            for (j=0; j<numCoords; j++)
                newClusters[index][j] += objects[i][j];
        }

        /* average the sum and replace old cluster center with newClusters */
        for (i=0; i<numClusters; i++) {
            for (j=0; j<numCoords; j++) {
                if (newClusterSize[i] > 0)
                    clusters[i][j] = newClusters[i][j] / newClusterSize[i];
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }

        delta /= numObjs;
    } while (delta > threshold && loop++ < 500);

    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);

    return 1;
}

void center_centroids_idx(int clust_rows, int clust_cols, float **cluster_centroids,
		int obj_rows, int obj_cols, float **objects, int *idx)
{

	float dist_frm_centroid[clust_rows];

	int i=0, j=0, l=0, m=0, b=0;

	for(i=0; i < obj_rows; i++){
		for(j=0; j < obj_cols; j++){

			for(l=0; l< clust_rows; l++){
				dist_frm_centroid[l] = 0;

				// find accumulated distance for each centroid
				for(m=0; m< clust_cols; m++){

					dist_frm_centroid[l] += sqrt((double)pow((objects[i][j] - cluster_centroids[l][m]), 2));

				}

			}
			printf("\nDIST ARRAY IN BEFORE..\n");
			for(b=0; b<clust_rows; b++)
				printf("\t %f", dist_frm_centroid[b]);

			idx[j+i*obj_cols] = min(dist_frm_centroid, clust_rows);

		}
	}

}

void center_centroids_idx2(int clust_rows, int clust_cols, float **cluster_centroids,
		int obj_rows, int obj_cols, float **objects, int *idx)
{

	float dist_frm_centroid[clust_rows];

	int i=0, j=0, l=0, m=0, b=0;

	for(i=0; i < obj_rows; i++){
		//for(j=0; j < obj_cols; j++){

			for(l=0; l< clust_rows; l++){
				dist_frm_centroid[l] = 0;

				// find accumulated distance for each centroid
				for(m=0; m< clust_cols; m++){

					dist_frm_centroid[l] += (double)pow((objects[i][m] - cluster_centroids[l][m]), 2);

				}

				dist_frm_centroid[l] = sqrt(dist_frm_centroid[l]);

			}

//			printf("\nDIST ARRAY IN BEFORE..\n");
//			for(b=0; b<clust_rows; b++)
//				printf("\t %f", dist_frm_centroid[b]);

			idx[i] = min(dist_frm_centroid, clust_rows);

		//}
	}

}

int min(float arr[], int length) {
    // returns the minimum value of array
	int i, idx = 0, b;

    double minimum = arr[0];

//    printf("\nDIST ARRAY IN FUNC..\n");
//    for(b=0; b<length; b++)
//    				printf("\t %f", arr[b]);

    for (i = 1; i < length; ++i) {
    	double temp = arr[i];
        if (minimum > arr[i]) {
            minimum = arr[i];
            idx = i;
        }
    }

    return idx;
}
