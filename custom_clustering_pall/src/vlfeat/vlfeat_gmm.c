/*
 ============================================================================
 Name        : Gmm_Test.c
 Author      : Vinay Gavirangaswamy
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <vlfeat/gmm.h>
#include <vlfeat/host.h>
#include <vlfeat/kmeans.h>
#include <vlfeat/vlad.h>
#include <vlfeat/generic.h>

#define TYPE doulbe
#define VL_F_TYPE VL_TYPE_DOUBLE



int max(double arr[], int length) {
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

void matrix_transpose(double *matrix, double *transpose, int row, int col){
	int c, d;
	   for (c = 0; c < row; c++)
	      for( d = 0 ; d < col ; d++ )
	         transpose[d * row + c] = matrix[c * col + d];
}


void vlfeat_gmm_adapter(int nclusters, int nrows,
		int ncols, double* data, int clusterid[])
/* Perform kmedoids clustering on data */
{

	int i,j;
	double * means;
	double * covariances;
	double * priors;
	double * posteriors, *transpose;
	double loglikelihood;

	vl_size numData = nrows;
	vl_size dimension = ncols;
	vl_size numClusters = nclusters;
	//vl_size maxiter = pmaxiter;
	//vl_size maxrep = pmaxrep;
	vl_size dataIdx, d, cIdx;



	VlRand rand;
	VlGMM * gmm;

//	printf("nrows = %d ncols=%d", nrows, ncols);
//
//	for(i=0; i<nrows; i++){
//		for(j=0;j<ncols; j++){
//
//			printf("%f ",data[j+i*(int)ncols]);
//
//		}
//		printf("\n");
//	}

	// create a new instance of a GMM object for float data
	gmm = vl_gmm_new(VL_TYPE_FLOAT, dimension, numClusters);


	vl_gmm_cluster (gmm, data, numData);
	posteriors = vl_gmm_get_posteriors(gmm);

	transpose = vl_malloc(sizeof(double) * numData * numClusters);

//	printf("\nPosteriors:\n");
//
//	for(i=0;i<numClusters; i++){
//		for(j=0; j< numData; j++)
//			printf("%f ", ((double*)posteriors)[i*numData+j]);
//		printf("\n");
//	}

	matrix_transpose(posteriors, transpose, numClusters, numData);

//	printf("\nTranspose:\n");
//
//	for(i=0;i<numData; i++){
//			for(j=0; j< numClusters; j++)
//				printf("%f ", ((double*)transpose)[i*numClusters+j]);
//			printf("\n");
//		}

	for(i=0;i<numData; i++){

		float temp [numClusters];

		memcpy(temp, &(transpose[i*numClusters]), sizeof(double) * numClusters);
		clusterid[i] = max(temp,numClusters);

//		printf("%d %d\n", max(temp,numClusters),clusterid[i]);

	}

//	printf("\nHELLO FROM GMM\n");
}
