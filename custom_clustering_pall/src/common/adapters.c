/*
 ============================================================================
 Name        : adapters.c
 Author      : Vinay B Gavirangaswamy
 Created on	 : Jul 26, 2015
 Version     : 1.0
 Copyright   :  This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.


    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.
 Description : 
 ============================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h> /* needed to define getpid() */
#include <sys/types.h>
#include <string.h>
#include <sys/errno.h>
#include <fcntl.h>
#include "common/wrapper.h"
#include "common/wrapperFunc.h"
#include "univ_of_tokyo/cluster.h"



/*
 * Adapter return an array of cluster centers of size [numClusters][numCoords] from library written by
 * Michiel de Hoon
 * University of Tokyo, Institute of Medical Science
 * Human Genome Center, Laboratory of DNA Information Analysis
 * Currently at RIKEN Genomic Sciences Center
 * mdehoon 'AT' gsc.riken.jp
 *
 */
int kmeans_adapter(int nclusters, int nrow, int ncol, double* data1d, int npass,
		char method, char dist, int clusterid[]) {
	int** mask;
	double* weight = malloc(ncol * sizeof(double));
	const int transpose = 0;
	int ifound = 0;
	double error;
	int i, j;
	double **data2d;

	mask = (int**) malloc(nrow * sizeof(int*));
	assert(mask != NULL);
	mask[0] = (int*) malloc(nrow * ncol * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrow; i++)
		mask[i] = mask[i - 1] + ncol;

	for (i = 0; i < nrow; i++)
		for (j = 0; j < ncol; j++)
			mask[i][j] = 1;

	/* allocate a 2D space for data as x[] (coordinates of cluster centers)
	 this array should be the same across all processes                  */
	data2d = (double**) malloc(nrow * sizeof(double*));
	data2d[0] = (double*) malloc(nrow * ncol * sizeof(double));
	for (i = 1; i < nrow; i++)
		data2d[i] = data2d[i - 1] + (ncol);

	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {

			data2d[i][j] = data1d[j + i * (int) ncol];

		}
	}

	kcluster(nclusters, nrow, ncol, data2d, mask, weight, transpose, npass,
			method, dist, clusterid, &error, &ifound);
}

int kclust_adapter2D(int nclusters, int nrow, int ncol, double** data2d,
		int npass, char* _method, char* _dist, int clusterid[]) {
	int** mask;
	double* weight = malloc(ncol * sizeof(double));
	const int transpose = 0;
	int ifound = 0;
	double error;
	int i, j;
	char method, dist;

	if (strcmp(_method, CNTR_FUN_MEAN) == 0) {
		method = 'a';
	} else if (strcmp(_method, CNTR_FUN_MEDIAN) == 0) {
		method = 'm';
	}

	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	} else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	} else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	} else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	} else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	mask = (int**) malloc(nrow * sizeof(int*));
	assert(mask != NULL);
	mask[0] = (int*) malloc(nrow * ncol * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrow; i++)
		mask[i] = mask[i - 1] + ncol;

	for (i = 0; i < nrow; i++)
		for (j = 0; j < ncol; j++)
			mask[i][j] = 1;

	for(j=0; j< ncol; j++){
		weight[j] = 1;
	}

//	printf("\nCalling kcluster(%d,%d, %d, %d, %d, %c, %c,%d, %d )...\n",nclusters, nrow,  ncol, transpose, npass, method, dist, error, ifound);

//	for(i=0; i< nrow; i++){
//		for(j=0; j< ncol; j++){
//			printf("%f ", data2d[i][j]);
//		}
//		printf("\n");
//	}

	kcluster(nclusters, nrow, ncol, data2d, mask, weight, transpose, npass,
			method, dist, clusterid, &error, &ifound);

	free(weight);

	free(mask[0]);
	free(mask);

}

/* ========================================================================= */

void agglom_adapter1D(int nclusters, int nrows, int ncols, double* data1d,
		char* _method, char* _dist, int *clusterid)
/* Perform hierarchical clustering on data */
{

	char method, dist;
	int i , j;


	// Assign lib specific parameters to link method
	// and distance function
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	} else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	} else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	} else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	} else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	} else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'c'; // TODO: ward linkage should be implemented later
	} else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	} else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'e';
	} else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	} else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	} else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	} else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	} else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	} else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	} else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	} else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}

	agglom_adapter(nclusters>= nrows ? nrows: nclusters,nrows,ncols,data1d,method,dist,clusterid);


}

void agglom_adapter2D(int nclusters, int nrows, int ncols, double** data2d,
		char* _method, char* _dist, int *clusterid)
/* Perform hierarchical clustering on data */
{

	char method, dist;



	// Assign lib specific parameters to link method
	// and distance function
	if (strcmp(_method, LNK_CODE_AVG) == 0) {
		method = 'a';
	} else if (strcmp(_method, LNK_CODE_CEN) == 0) {
		method = 'c';
	} else if (strcmp(_method, LNK_CODE_COM) == 0) {
		method = 'm';
	} else if (strcmp(_method, LNK_CODE_SIN) == 0) {
		method = 's';
	} else if (strcmp(_method, LNK_CODE_MED) == 0) {
		method = 'a'; // TODO: median linkage should be implemented later
	} else if (strcmp(_method, LNK_CODE_WAR) == 0) {
		method = 'c'; // TODO: ward linkage should be implemented later
	} else if (strcmp(_method, LNK_CODE_WEI) == 0) {
		method = 'a'; // TODO: weighted linkage should be implemented later
	}

	if (strcmp(_dist, DIST_MTRC_EUC) == 0) {
		dist = 'e';
	} else if (strcmp(_dist, DIST_MTRC_SEU) == 0) {
		dist = 'e';
	} else if (strcmp(_dist, DIST_MTRC_CIT) == 0) {
		dist = 'b';
	} else if (strcmp(_dist, DIST_MTRC_COR) == 0) {
		dist = 'c';
	} else if (strcmp(_dist, DIST_MTRC_COS) == 0) {
		dist = 'o';
	} else if (strcmp(_dist, DIST_MTRC_MAH) == 0) {
		dist = 'm';
	} else if (strcmp(_dist, DIST_MTRC_JAC) == 0) {
		dist = 'j';
	} else if (strcmp(_dist, DIST_MTRC_CHE) == 0) {
		dist = 'h';
	} else if (strcmp(_dist, DIST_MTRC_SPE) == 0) {
		dist = 's';
	} else if (strcmp(_dist, DIST_MTRC_HAM) == 0) {
		dist = 'g';
	}


	agglom_adapter_2d(nclusters>= nrows ? nrows: nclusters,nrows,ncols,data2d,method,dist,clusterid);


}

void agglom_adapter_2d(int nclusters, int nrows, int ncols, double **data2d,
		char method, char dist, int *clusterid)
/* Perform hierarchical clustering on data */
{
	int i, j;

	int** mask;

	double* weight = malloc(ncols * sizeof(double));

	//const int nnodes = nrows - 1;
	Node* tree;

	for (i = 0; i < ncols; i++)
		weight[i] = 1.0;

	mask = (int**) malloc(nrows * sizeof(int*));
	assert(mask != NULL);
	mask[0] = (int*) malloc(nrows * ncols * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrows; i++)
		mask[i] = mask[i - 1] + ncols;



	for (i = 0; i < nrows; i++)
				for (j = 0; j < ncols; j++)
					if (data2d[i][j] != 0)
						mask[i][j] = 1;
					else
						mask[i][j] = 0;


	tree = treecluster(nrows, ncols, data2d, mask, weight, 0, dist, method, 0);
	if (!tree) { /* Indication that the treecluster routine failed */
		printf("treecluster routine failed due to insufficient memory\n");
		free(weight);
		return;
	}

	/*printf("Node     Item 1   Item 2    Distance\n");
	 for(i=0; i<nnodes; i++)
	 printf("%3d:%9d%9d      %g\n",
	 -i-1, tree[i].left, tree[i].right, tree[i].distance);
	 printf("\n");*/
	//printf("=============== Cutting a hierarchical clustering tree ==========\n");
	//clusterid = malloc(nrows*sizeof(int));

	cuttree(nrows, tree, nclusters, clusterid);
//	for(i=0; i<nrows; i++)
//	 printf("obj[%2d]: cluster %2d\n", i, clusterid[i]);
//	 printf("\n");

	free(tree);
	free(data2d[0]);
	free(data2d);
	free(mask[0]);
	free(mask);
	free(weight);

	return;
}

void agglom_adapter(int nclusters, int nrows, int ncols, double* data1d,
		char method, char dist, int *clusterid)
/* Perform hierarchical clustering on data */
{
	int i, j;

	int** mask;
	double **data2d;
	double* weight = malloc(ncols * sizeof(double));

	//const int nnodes = nrows - 1;
	Node* tree;

	for (i = 0; i < ncols; i++)
		weight[i] = 1.0;

	mask = (int**) malloc(nrows * sizeof(int*));
	assert(mask != NULL);
	mask[0] = (int*) malloc(nrows * ncols * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrows; i++)
		mask[i] = mask[i - 1] + ncols;


	/* allocate a 2D space for data as x[] (coordinates of cluster centers)
	 this array should be the same across all processes                  */
	data2d = (double**) malloc(nrows * sizeof(double*));
	data2d[0] = (double*) malloc(nrows * ncols * sizeof(double));
	for (i = 1; i < nrows; i++)
		data2d[i] = data2d[i - 1] + (ncols);

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			data2d[i][j] = data1d[i*ncols+j];
		}
	}

	for (i = 0; i < nrows; i++)
				for (j = 0; j < ncols; j++)
					if (data2d[i][j] != 0)
						mask[i][j] = 1;
					else
						mask[i][j] = 0;


	tree = treecluster(nrows, ncols, data2d, mask, weight, 0, dist, method, 0);
	if (!tree) { /* Indication that the treecluster routine failed */
		printf("treecluster routine failed due to insufficient memory\n");
		free(weight);
		return;
	}

	/*printf("Node     Item 1   Item 2    Distance\n");
	 for(i=0; i<nnodes; i++)
	 printf("%3d:%9d%9d      %g\n",
	 -i-1, tree[i].left, tree[i].right, tree[i].distance);
	 printf("\n");*/
	//printf("=============== Cutting a hierarchical clustering tree ==========\n");
	//clusterid = malloc(nrows*sizeof(int));

	cuttree(nrows, tree, nclusters, clusterid);
//	for(i=0; i<nrows; i++)
//	 printf("obj[%2d]: cluster %2d\n", i, clusterid[i]);
//	 printf("\n");

	free(tree);
	free(data2d[0]);
	free(data2d);
	free(mask[0]);
	free(mask);
	free(weight);

	return;
}

/* ========================================================================= */

void kmedoids_adapter(int nclusters, int npass, int nrows, int ncols,
		double* data1d, char method, char dist, int clusterid[])
/* Perform kmedoids clustering on data */
{
	int ifound = 0;
	double error;
	int i, j;
	double** distmatrix;
	int** mask;
	double **data2d;
	double* weight = malloc(ncols * sizeof(double));

//	  printf("\nValue of ncols %d\n", ncols);

	for (i = 0; i < ncols; i++)
		weight[i] = 1.0;
//
	mask = (int**) malloc(nrows * sizeof(int*));
	assert(mask != NULL);

	mask[0] = (int*) malloc(nrows * ncols * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrows; i++)
		mask[i] = mask[i - 1] + ncols;

	for (i = 0; i < nrows; i++)
		for (j = 0; j < ncols; j++)
			mask[i][j] = 1;

	/* allocate a 2D space for data as x[] (coordinates of cluster centers)
	 this array should be the same across all processes                  */
	data2d = (double**) malloc(nrows * sizeof(double*));
	data2d[0] = (double*) malloc(nrows * ncols * sizeof(double));
	for (i = 1; i < nrows; i++)
		data2d[i] = data2d[i - 1] + (ncols);

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {

			data2d[i][j] = data1d[j + i * (int) ncols];

		}
	}

	distmatrix = distancematrix(nrows, ncols, data2d, mask, weight, dist, 0);

	kmedoids(nclusters, nrows, distmatrix, npass, clusterid, &error, &ifound);

}



/*
 * Utility function that used dist function from university of tokyo library
 */

double** getDistMetric (int nrows, int ncols,
		double** data, char dist){

	int i, j;
	int** mask;

	double* weight = malloc(ncols * sizeof(double));

	for (i = 0; i < ncols; i++)
		weight[i] = 1.0;
//
	mask = (int**) malloc(nrows * sizeof(int*));
	assert(mask != NULL);

	mask[0] = (int*) malloc(nrows * ncols * sizeof(int));

	assert(mask[0] != NULL);
	for (i = 1; i < nrows; i++)
		mask[i] = mask[i - 1] + ncols;

	for (i = 0; i < nrows; i++)
		for (j = 0; j < ncols; j++)
			if (data[i][j] != 0)
				mask[i][j] = 1;
			else
				mask[i][j] = 0;

	  /* Set the metric function as indicated by dist */
	return distancematrix(nrows, ncols, data, mask, weight, dist, 0);

}
