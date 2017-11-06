/*
 ============================================================================
 Name        : spectralWrapper.c
 Author      : Vinay B Gavirangaswamy
 Created on	 : Aug 28, 2015
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


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>     /* strtok() */
#include <cerrno>
#include <cassert>
#include <ctime>
#include <unistd.h>
#include <fstream>
#include <mpi.h>


extern "C" {
#include "common/wrapperFunc.h"
void file_read_dim(
                  char *filename,      /* input file name */
                  int  *numObjs,       /* no. data objects (local) */
                  int  *numCoords);     /* no. coordinates */
}
#include "common/wrapper.h"
#include <common/pSpectral.h>
#include "pspectral/compute_distance.h"
#include "pspectral/distance_to_similarity.h"
#include "pspectral/evd.h"
#include "pspectral/kmeans.h"


/* =========================================================================
 *
 * Spectral clustering wrapper
 *
 * ========================================================================= */

void spectral_adapter2D(MPI_Comm new_comm, int nclusters, int nrows, int ncols, double** data2d,
		char* _method, char* _dist, int *clusterid/*, char* tempDataInputFile,
		char* tempDistOutputFile, char* tempDist2SimFile, char* tempEvdValsFile,
		char* tempEvdVecFile, char* tempClustIdxFile*/){

	char dist;
	int myrank, nproc;


	MPI_Comm_size(new_comm, &nproc);
	MPI_Comm_rank(new_comm, &myrank);



		// Assign lib specific parameters to link method
		// and distance function
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
		} else{
			dist = 'e';
		}



	internalOrchestrator(new_comm, data2d, nrows, ncols, dist, clusterid, nclusters/*, tempDataInputFile,
			 tempDistOutputFile, tempDist2SimFile, tempEvdValsFile, tempEvdVecFile, tempClustIdxFile*/);
	d_printf("Returned from internalOrchestrator...[My rank is %d]\n", myrank);

}



/* =========================================================================
 *
 * Spectral clustering wrapper
 *
 * ========================================================================= */
using namespace learning_psc;
using namespace std;

int syncAllProcess(int myrank, int& token, MPI_Comm new_comm) {
	if (myrank == 0) {
		token = std::rand() % 99 + 1;

		MPI_Bcast(&token, 1, MPI_INT, 0, new_comm);

	} else {

//		d_printf("[MY-RANK %d TOKEN %d] Waiting on master process \n", myrank,	token);
		MPI_Bcast(&token, 1, MPI_INT, 0, new_comm);

	}
	return myrank;
}

/* ========================================================================= */
void internalOrchestrator(MPI_Comm new_comm, double **data2d, int nrows, int ncols, char dist,
		int *clusterid, int nclusters /*,char* tempDataInputFile,
		char* tempDistOutputFile, char* tempDist2SimFile, char* tempEvdValsFile,
		char* tempEvdVecFile, char* tempClustIdxFile*/) {

	int myrank, nproc;
	int token;
	int evdVec_count;

	// variable for debugging
	char *line;
	size_t len;
	ssize_t numBytesRead;
	FILE *ifp;

	std::stringstream tempDataInputBuffer;

	const char *dataInputPtr, *datComputeDistancePtr, *dataCompute_d2sPtr;

	MPI_Comm_size(new_comm, &nproc);
	MPI_Comm_rank(new_comm, &myrank);

	if (myrank == 0) {

//		std::ofstream ofp;
//		ofp .open(tempDataInputFile);

			for(int i=0;i<nrows; i++){
				for(int j=0;j<ncols; j++){
//					ofp<<j<<":"<< data2d[i][j] << " ";
					tempDataInputBuffer<<j<<":"<< data2d[i][j] << " ";
//					printf("%d:%2.2f ", j, data2d[i][j]);
				}
//				ofp<<std::endl;
				tempDataInputBuffer<<std::endl;
//				printf("\n");
			}

//			ofp.close();


			token = std::rand() % 99 + 1;

			MPI_Bcast(&token, 1, MPI_INT, 0, new_comm);

		} else {

//			d_printf("[MY-RANK %d TOKEN %d] Waiting on master process \n",myrank, token);
			MPI_Bcast(&token, 1, MPI_INT, 0, new_comm);

		}

//	MPI_Barrier(MPI_COMM_WORLD);

//	d_printf("[MY-RANK %d TOKEN %d] Starting Step 1: Compute Similarity  dist = %c \n", myrank, token, dist);



	/*Step 1: Compute Similarity*/

	ComputeDistance compute(10, new_comm);

	dataInputPtr = new char[tempDataInputBuffer.str().size()+1];
	strcpy(dataInputPtr, tempDataInputBuffer.str().c_str());
	compute.Read(dataInputPtr);

//	d_printf("[MY-RANK %d TOKEN %d] Using euclidean distance dist = %c \n", myrank, token, dist);

	  if (dist == 'e') {
			// do euclidean distance
		  compute.ParallelComputeEuclidean();
		} else if (dist== 'b') {
			// do city block distance
			compute.ParallelComputeManhattan();

		} else if (dist== 'c') {
			// do correlation distance
			compute.ParallelComputeCorrelation();

		} else if (dist== 'o') {
			// do cosine distance
			compute.ParallelComputeCosine();

		}else {
			// do euclidean distance by default
			compute.ParallelComputeEuclidean();
		}


//	  datComputeDistancePtr = compute.Write(tempDistOutputFile);
	  datComputeDistancePtr = compute.Write();


//	  d_printf(" RETURNED FOMR COMPUTE DISTANCE STEP:\n %s \n", datComputeDistancePtr);

	  myrank = syncAllProcess(myrank, token, new_comm);

	  /*Step 2: Compute distance_to_similarity*/


	  DistanceToSimilarity compute_d2s(new_comm);
//	  compute_d2s.ReadAndSetSymmetric(tempDistOutputFile);
	  compute_d2s.ReadAndSetSymmetric(datComputeDistancePtr);
	  compute_d2s.Compute();
//	  compute_d2s.Write(tempDist2SimFile);
	  dataCompute_d2sPtr = compute_d2s.Write();



	  myrank = syncAllProcess(myrank, token, new_comm);

//	  d_printf(" RETURNED FOMR DISTANCE TO SIMILARITY  STEP:\n %s \n", dataCompute_d2sPtr);


	  /*Step 3: Compute evd*/

	  int FLAGS_eigenvalue = 100;
	  int FLAGS_eigenspace = 300;
	  int FLAGS_arpack_iterations = 300;
	  double FLAGS_arpack_tolerance = 0.1;
//	  string FLAGS_input(tempDist2SimFile);
//	  string FLAGS_eigenvalues_output(tempEvdValsFile);
//	  string FLAGS_eigenvectors_output(tempEvdVecFile);

	  char /**dataEigenvalues_output,*/ *eigenvectors_output;

	  Evd evd(new_comm);
//	  evd.Read(FLAGS_input);
	  evd.Read(dataCompute_d2sPtr);
	  evd.Compute(FLAGS_eigenvalue, FLAGS_eigenspace, FLAGS_arpack_iterations, FLAGS_arpack_tolerance);


//	  evd.Write(FLAGS_eigenvalues_output, FLAGS_eigenvectors_output);
	  eigenvectors_output = evd.Write(/*dataEigenvalues_output,*/ );


	myrank = syncAllProcess(myrank, token, new_comm);

//	d_printf(" RETURNED FROM EVD  STEP:\n %s \n", eigenvectors_output);

//	d_printf("FLAGS_eigenvectors_output length:%d \n", evdVec_count);

	  /*Step 4: Compute kmeans on evd for idx*/


	  int FLAGS_num_clusters = nclusters;
	  int FLAGS_kmeans_loop  = 100;
	  double FLAGS_kmeans_threshold = 1e-3;
	  string FLAGS_initialization_method = "orthogonal_centers";
	  char *datakmeans_output;
//	  string FLAGS_kmeans_output(tempClustIdxFile);

	  KMeans kmeans(new_comm);
//	  kmeans.Read(FLAGS_eigenvectors_output);
	  kmeans.Read(eigenvectors_output);
	  kmeans.DoKMeans(FLAGS_num_clusters, FLAGS_initialization_method,
	                    FLAGS_kmeans_loop, FLAGS_kmeans_threshold);
//	  kmeans.Write(FLAGS_kmeans_output);
	  datakmeans_output = kmeans.Write();

//	  d_printf(" RETURNED FROM KMEANS  STEP:\n %s \n", datakmeans_output);

	  myrank = syncAllProcess(myrank, token, new_comm);

	if (myrank == 0) {

		std::istringstream fin2(datakmeans_output);
		 string line;
		  int i = 0;
		  while (getline(fin2, line) && i < nrows) {  // Each line is a training document.
		    if (line.size() > 0 &&      // Skip empty lines.
		        line[0] != '\r' &&      // Skip empty lines.
		        line[0] != '\n' &&      // Skip empty lines.
		        line[0] != '#') {       // Skip comment lines.

		    	clusterid[i] = atoi(line.c_str());
//		    	std::cout<<" "<< clusterid[i];
		    	i++;
		    }
		  }

//		FILE * ifp;
//		ifp = fopen(tempClustIdxFile, "r");
//
//		for (int i = 0; i < nrows; i++) {
//			if (!fscanf(ifp, "%d", &clusterid[i]))
//				break;
////					       std::cout<<" "<< clusterid[i];
//		}
////		 std::cout<<std::endl;
//
//		fclose(ifp);
	} else {
		for (int i = 0; i < nrows; i++) {
			clusterid[i] = NaN;
			//					       std::cout<<" "<< clusterid[i];
		}

	}

	free(dataInputPtr);
	free(datakmeans_output);
	free(eigenvectors_output);
	free(dataCompute_d2sPtr);
	free(datComputeDistancePtr);


//	if (myrank == 0) {
//			token = std::rand() % 99 + 1;
//			for(int i=0; i< nproc; i++)
//		    MPI_Send(&token, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
//		} else if (myrank != 0) {
//
//		    MPI_Recv(&token, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
//		             MPI_STATUS_IGNORE);
//	}

	MPI_Barrier(new_comm);

}
