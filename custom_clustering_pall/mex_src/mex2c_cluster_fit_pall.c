/*
 * mex2c_cluster_fit.c
 * CLUSTER_FIT - Clusters data using specified method and parameters
 *
 * Syntax:  [C] = cluster_fitmodel(X,bootIDxs,kList,models,PST)
 * INPUTS ARGUMENTS
 *			  X  :  [n,p] data (after preprocessing, if preprocessed)
 *		bootIDxs :  [N,B] array of bootstrap indices for data X
 *					specify [] for no boostrapping.
 *		   kList :  array of values of k to be fit; can be scalar
 * 	 	  models :  clustering method to be fit e.g. {'agg','spe','med','kme'}
 *	 	 	 PST :  parameter structure
 *
 * VARIABLE ARGUMENTS
 *
 * OUTPUT ARGUMENTS
 * 			 C   :  clustering result output structure
 *
 *  Created on: Aug 10, 2015
 *      Author: Vinay B Gavirangaswamy
 *  Copyright   :  This program is free software: you can redistribute it and/or modify
 * 	    		   it under the terms of the GNU General Public License as published by
 *		   		   the Free Software Foundation, either version 3 of the License, or
 *   			   (at your option) any later version.
 *
 *   			   This program is distributed in the hope that it will be useful,
 *   			   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   			   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   			   GNU General Public License for more details.
 *
 *
 *   			   You should have received a copy of the GNU General Public License
 *   			   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "mex.h"
#include <stdlib.h> /* needed to define exit() */
#include <unistd.h> /* needed to define getpid() */
#include <stdio.h> /* needed for printf() */
#include <sys/types.h>
#include <assert.h>
#include <sys/errno.h>
#include <fcntl.h>
#include "common/wrapper.h"





void doClusterFitSeq(char *modelListFile, char *kListFile, char *bootIdxFilesFile,
		char *final_output_file,
		int nBootstrap, int nRow, int nCol, int nModels,
		int num_mpi_threads) {

	pid_t childpid;
	char num_t[32], s_nBootstrap[32], num_features[32], s_nModels[32], s_nRows[32];
	int count = 0;

	char *config_str, *executalbe;

	sprintf(num_t, "%d", num_mpi_threads);
	sprintf(s_nBootstrap, "%d", nBootstrap);
	sprintf(num_features, "%d", nCol);
	sprintf(s_nModels, "%d", nModels);
	sprintf(s_nRows, "%d", nRow);

	executalbe = malloc(TMP_FILE_NM_SIZE * sizeof(char));

	strcpy(executalbe, SEQ_CLUST_ALGO_EXECUTABLE);



	/*Spawn a child to run the program.*/
	childpid = fork();
	if (childpid == 0) { /* child process */
		/*char *av[] = { mpiReadFd, mpiWriteFd, (char *) 0 }; /* each element represents a command line argument */

		/*execv("Matlab_MPI", av);*/

		char **args = (char **) calloc(19, sizeof(char *));
		args[0] = "mpiexec";
		args[1] = "-n";
		args[2] = num_t;
		args[3] = executalbe;"../seqClustAlgo";
		args[4] = "-m";
		args[5] = modelListFile;
		args[6] = "-n";
		args[7] = s_nModels;
		args[8] = "-o";
		args[9] = final_output_file;
		args[10] = "-c";
		args[11] = num_features;
		args[12] = "-b";
		args[13] = bootIdxFilesFile;
		args[14] = "-s";
		args[15] = s_nBootstrap;
		args[16] = "-r";
		args[17] = s_nRows;
		args[18] = (char *) 0;

		d_printf("MPI argument...\n");

		for(count = 0; count < 19; count++)
			printf("%s ", args[count]);

		printf("\n");

		/*args[10] = (char *) 0;*/

		execvp("mpirun", args);
		perror("execvp"); /* if we get here, execvp failed */

	} else { /* pid!=0; parent process */
		/*printf("\nParent continues..\n");*/
		waitpid(childpid, 0, 0); /* wait for child to exit */
	}
}


/**
 * @Desction:
 *
 * @Input: X_IN should be of type double data type
 * 		   PST_IN should be of type structure
 * 		   BOOTIDX_IN should be of type matlab uint16
 */

/* Always include this */
void mexFunction(int nlhs, mxArray *plhs[],
/* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{

#define X_IN			prhs[0]
#define PST_IN			prhs[1]
#define BOOTIDX_IN		prhs[2]
#define FINAL_OUTPUT_FILE prhs[3]
#define NUM_THREADS		prhs[4]


	int num_mpi_threads = 0;
	mxArray *mxCellArrPtr;
	int ifield, nfields, field_num;
	mxClassID *classIDflags;
	mwSize nStructElems;
	mxArray *tmp, *pst_ensemble;
	const mxArray *cell_element_ptr;
	char *name;
	int VERBOSE_Q, INCLUDE_REPsQ, INCLUDE_CENTERsQ;
	double *kList;
	char **clustList;
	char **modelList;
	int nBootstrap;
	int nReps;
	int n_klist, n_modellist;
	int kmeanModelCnt = 0, aggloModelCnt = 0, gmmModelCnt = 0, specModelCnt = 0;
	int c_nModels, c_nKList, c_nBootstrap, c_nReps;
	/*
	 * variables used to make function call to matlab functions
	 */
	double *x_actual_data, *x_transpose;
	int m_x, n_x;
	int count_x_transpose;
	/*mxArray *_rhs[3], *_lhs[1], *tmp_rhs[3];
	 int _nlhs, _nrhs;*/
	int i, j, m, n;
	/*
	 * variable used for handling bootstrap data to file
	 */
	char **fnames; /* pointers to temp file names */
	char *modelListFile, *kListFile, *bootIdxFilesFile, *kmeanOutIdxFile,
			*gmmOutIdxFile, *aggloOutIdxFile, *specOutIdxFile;
	FILE *ifp, *ofp;
	unsigned short *bootIdxs;
	double *data_bootstrap;
	int m_data_bootstrap, n_data_bootstrap, count_m_data_bootstrap;
	double *bootIdxs_col;
	int m_bootIdxs, n_bootIdxs;

	int count = 0;

	int op_count;
	/*
	 * Variable to prepare output to matlab from function
	 *
	 */
	char *final_output_file;

	FILE *final_opfile;

	const char *field_names[] = { "partitions", "centers" };
	const char *partitions_field_names[] = { "model", "boot", "k", "rep",
			"indices" };
	const char *centers_field_names[] = { "model", "boot", "k", "rep",
			"centernum", "center" };

	int partitions_field, centers_field;

	char **mode_string;
	double ensembleRslt_rows = 0;



	/**
	 * Extract data from PST struct
	 */
	/* get input arguments */
	final_output_file = mxArrayToString(FINAL_OUTPUT_FILE);

	nfields = mxGetNumberOfFields(PST_IN);
	nStructElems = mxGetNumberOfElements(PST_IN);
	/* check proper input and output */
	if (nStructElems != 1)
		mexErrMsgIdAndTxt("MATLAB: Multiple PST structures",
				"Only one input is required.");


	tmp = mxGetField(PST_IN, 0, VERBOSE);
	VERBOSE_Q = mxGetScalar(tmp);
	field_num = mxGetFieldNumber(PST_IN, NAME);
	tmp = mxGetFieldByNumber(PST_IN, 0, field_num);
	name = mxArrayToString(tmp);

	pst_ensemble = mxGetField(PST_IN, 0, ENSEMBLE);

	num_mpi_threads = mxGetScalar(NUM_THREADS);

	/**
	 * Extract data from ENSEMBLE struct
	 * ENSEMBLE is inside PST
	 */
	INCLUDE_REPsQ = mxGetScalar(
			mxGetField(pst_ensemble, 0, ENSEMBLE_INCLUDEREPSQ));
	INCLUDE_CENTERsQ = mxGetScalar(
			mxGetField(pst_ensemble, 0, ENSEMBLE_INCLUDECENTERSQ));
	nBootstrap = mxGetScalar(mxGetField(pst_ensemble, 0, ENSEMBLE_NBOOTSTRAPS));

	if (INCLUDE_REPsQ == 1) {
		nReps = mxGetScalar(mxGetField(pst_ensemble, 0, ENSEMBLE_NREPS));
	} else {
		nReps = 1;
	}

	tmp = mxGetField(pst_ensemble, 0, ENSEMBLE_KLIST);
	n_klist = mxGetN(tmp);
	kList = mxGetPr(tmp);

/*	kListFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(kListFile);

	ofp = fopen(kListFile, "w");
	for (j = 0; j < n_klist; j++) {
		fprintf(ofp, "%d\n", (int) kList[j]);
	}
	fclose(ofp);*/


	/*
	 * Get list of model and cluster to run.
	 */
	tmp = mxGetField(pst_ensemble, 0, ENSEMBLE_MODELLIST);
	n_modellist = mxGetNumberOfElements(tmp);

	modelList = (char**) malloc(n_modellist * sizeof(char*));

	modelListFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(modelListFile);



	for (j = 0; j < n_modellist; j++) {
		cell_element_ptr = mxGetCell(tmp, j);
		modelList[j] = mxArrayToString(cell_element_ptr);
		/*fprintf(ofp, "%s\n", modelList[j]);*/
	}

	/*fclose(ofp);*/

	/* allocate memory  for storing pointers */

	clustList = (char**) malloc(CLUST_IN_ENSEMBLE * sizeof(char*));
	for (j = 0; j < CLUST_IN_ENSEMBLE; j++) {
		clustList[j] = malloc(3 * sizeof(char));
		memset(clustList[j], 0, sizeof(clustList[j]));
	}

	getClusteringMethodsList(n_modellist, modelList, clustList);

	/**
	 * Create temporary files to transfer
	 * data for clustering.
	 */
	/* allocate memory  for storing pointers */
	fnames = (char**) malloc(nBootstrap * sizeof(char*));
	for (i = 0; i < nBootstrap; i++) {
		fnames[i] = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	}

	m_bootIdxs = mxGetM(BOOTIDX_IN);
	n_bootIdxs = mxGetN(BOOTIDX_IN);
	bootIdxs = mxGetPr(BOOTIDX_IN);

	m_x = mxGetM(X_IN);
	n_x = mxGetN(X_IN);
	x_actual_data = mxGetPr(X_IN);

	if (m_x != m_bootIdxs) {
		mexPrintf(
				"Error: Base data (X) and bootstrap indices (bootIdxs) used to generate bootstrap data sample do not match.\n");
		return;
	}

	/*
	 * Transpose X as data stored in single dimensional array
	 * is row major not column
	 */
	x_transpose = malloc(m_x * n_x * sizeof(double));
	count_x_transpose = 0;
	for (i = 0; i < m_x; i++) {
		for (j = 0; j < n_x; j++) {
			x_transpose[count_x_transpose] = x_actual_data[i + j * m_x];
			count_x_transpose++;
		}
	}

	m_data_bootstrap = m_bootIdxs;
	n_data_bootstrap = n_x;

	data_bootstrap = malloc(
			m_data_bootstrap * n_data_bootstrap * sizeof(double));



	/**
	 * Generate bootstrap data and write it
	 * to a temporary file to use by clustering
	 * algorithms.
	 * Create bootstrap data sets...
	 */
	for (i = 0; i < n_bootIdxs; i++) {

		fnames[i] = malloc(TMP_FILE_NM_SIZE * sizeof(char));
		createTempFile(fnames[i]);

		for (j = 0; j < m_bootIdxs; j++) {
			/*printf("[%d %d] ",(int) bootIdxs[j + i * m_bootIdxs], (int)bootIdxs[j + i * m_bootIdxs]*n_x);*/
			memcpy(
					&data_bootstrap[j * n_data_bootstrap],
					&x_transpose[((int) bootIdxs[j + i * m_bootIdxs] - 1) * n_x],
					n_x * sizeof(double));
		}
		/*printf("\n");*/
		ofp = fopen(fnames[i], "w");
		for (m = 0; m < m_data_bootstrap; m++) {
			for (n = 0; n < n_data_bootstrap; n++) {

				fprintf(ofp, "%f ", data_bootstrap[n + m * n_data_bootstrap]);
			}

			fprintf(ofp, " \n");
		}
		fclose(ofp);
	}

	bootIdxFilesFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));

	createTempFile(bootIdxFilesFile);

	ofp = fopen(bootIdxFilesFile, "w");
	for (i = 0; i < n_bootIdxs; i++) {
		fprintf(ofp, "%s\n", fnames[i]);
	}
	fclose(ofp);



	/*for (thread_iter = 0; thread_iter < NUM_CLUSTERS_IN_ENSEMBLE; thread_iter++) {
	 printf("In main: creating thread %ld\n", thread_iter);
	 parameters.tid = thread_iter;
	 parameters.clustList = clustList;
	 parameters.modelListFile = modelListFile;
	 parameters.kListFile = kListFile;
	 rc = pthread_create(&threads[thread_iter], NULL, doClusterFit, (void *) &parameters);
	 if (rc) {
	 printf("ERROR; return code from pthread_create() is %d\n", rc);
	 exit(-1);
	 }
	 }
	 */
	/*
	 * Create temprory files to get clustering idx output
	 */
/*	kmeanOutIdxFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(kmeanOutIdxFile);

	gmmOutIdxFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(gmmOutIdxFile);

	aggloOutIdxFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(aggloOutIdxFile);

	specOutIdxFile = malloc(TMP_FILE_NM_SIZE * sizeof(char));
	createTempFile(specOutIdxFile);*/

	/*Calculate how many models in each clustering type*/

/*	kmeanModelCnt = getClustMthdCnt(KMEANS_SHRT, modelList, n_modellist);
	aggloModelCnt = getClustMthdCnt(AGGLO_SHRT, modelList, n_modellist);
	gmmModelCnt = getClustMthdCnt(GMM_SHRT, modelList, n_modellist);
	specModelCnt = getClustMthdCnt(SPECTRAL_SHRT, modelList, n_modellist);*/



	/*Do the clustering*/
	mexPrintf("\n---------------------------------\n");
	mexPrintf("Generating cluster ensemble\n");
	mexPrintf("---------------------------------\n");
	mexPrintf("\nRunning Ensemble with %d Bootstrap...\n", nBootstrap);
	int total_models = 1;
	ofp = fopen(modelListFile, "w");
	for (c_nModels = 0; c_nModels < n_modellist; c_nModels++) {
		for (c_nKList = 0; c_nKList < n_klist; c_nKList++) {
			for (c_nBootstrap = 0; c_nBootstrap < nBootstrap; c_nBootstrap++) {
				for (c_nReps = 1; c_nReps < nReps+1; c_nReps++) {

					total_models++;
					/*printf("\n%s %.0f %s %d", modelList[c_nModels], kList[c_nKList], fnames[c_nBootstrap], c_nReps);*/
					fprintf(ofp,"%s %.0f %d %d\n", modelList[c_nModels], kList[c_nKList], c_nBootstrap, c_nReps);


				}
			}

		}
	}

	fclose(ofp);

	if (final_opfile = fopen(final_output_file, "r")) {
		fclose(final_opfile);
		remove(final_output_file);
	}


	doClusterFitSeq(modelListFile, kListFile, bootIdxFilesFile, final_output_file, nBootstrap,m_data_bootstrap, n_data_bootstrap,
			total_models, num_mpi_threads);
/*
	printf("\nRunning p-spectral now...\n");
	doClusterFitPSpec(modelListFile, kListFile, bootIdxFilesFile, specOutIdxFile, nBootstrap, nReps,
			n_modellist, n_klist, num_mpi_threads);
*/

	/*
	 * Read cluster idx data from corresponding output file
	 * and populate output matlab struct()
	 * kmeans: 		kmeanOutIdxFile
	 * gmm: 		gmmOutIdxFile
	 * agglo: 		aggloOutIdxFile
	 * spectral: 	specOutIdxFile
	 */

	/*ensembleRslt_rows = n_modellist * TMP_NKLIST * nBootstrap * TMP_NREPS;
	ensembleRslt_rows = ((RUN_KMEANS == 1 ? kmeanModelCnt : 0)
			+ (RUN_AGGLO == 1 ? aggloModelCnt : 0)+ (RUN_SPECTRAL == 1 ? specModelCnt : 0)) * n_klist * nBootstrap
			* nReps;

	op_count = 0;
*/	mexPrintf("FINAL_OUTPUT_FILE = %s", final_output_file);


	/*
	 * If PST has kmeans extract idx from file to return
	 */
/*	if (doKmeans(clustList) && RUN_KMEANS) {
		fetchNpoulateRslts(KMEANS_SHRT, clustList, modelList, n_modellist,
				kmeanOutIdxFile, nBootstrap, kList, bootIdxs, m_bootIdxs,
				final_opfile, n_klist, nReps);
	}

	 * If PST has Agglo extract idx from file to return

	if (doAgglom(clustList) && RUN_AGGLO) {

		fetchNpoulateRslts(AGGLO_SHRT, clustList, modelList, n_modellist,
				aggloOutIdxFile, nBootstrap, kList, bootIdxs, m_bootIdxs,
				final_opfile, n_klist, nReps);
	}

	 * If PST has Spectral extract idx from file to return

	if (doSpectral(clustList) && RUN_SPECTRAL) {

		fetchNpoulateRslts(SPECTRAL_SHRT, clustList, modelList, n_modellist,
				specOutIdxFile, nBootstrap, kList, bootIdxs, m_bootIdxs,
				final_opfile, n_klist, nReps);
	}*/

	/*Unlink unnecessary files*/

	for (i = 0; i < nBootstrap; i++) {
		unlink(fnames[i]);
	}

	unlink(modelListFile);

	unlink(bootIdxFilesFile);

//	unlink(kListFile);
//	unlink(kmeanOutIdxFile);
//	unlink(gmmOutIdxFile);
//	unlink(aggloOutIdxFile);
//	unlink(specOutIdxFile);


	mexPrintf("\n");

	/* Last thing that main() should do */
	/*pthread_exit(NULL);*/

}

