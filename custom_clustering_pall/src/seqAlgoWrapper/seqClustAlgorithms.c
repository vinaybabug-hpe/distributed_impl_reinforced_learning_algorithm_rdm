/*
 * run_clust_p_mchine.c
 * CLUSTER_FIT - Clusters data using specified method and parameters
 * 				 on 1 cluster per machine.
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
 *  Created on: Aug 21, 2015
 *      Author: Vinay B Gavirangaswamy
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
 *
 */
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* getopt() */
#include <mpi.h>
#include <errno.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <malloc.h>

extern "C" {
void file_read_dim(
                  char *filename,       /* input file name */
                  int  *numObjs,       /* no. data objects (local) */
                  int  *numCoords);     /* no. coordinates */
#include "common/wrapperFunc.h"
}
#include "common/wrapper.h"
#include <common/pSpectral.h>

int _debug;

/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0) {
	char *help =
			"Usage: %s [switches] -i filename -n num_clusters\n"
					"       -m modelListFile    : file containing models used on data to be clustered\n"
					"       -n centers     		: total number of models used on data to be clustered\n"
					"       -k kListFile   		: file containing list of K cluster to produce on data to be clustered\n"
					"       -i filename    		: filename containing names of file that has bootstrap data\n"
					"       -b #bootstraps		: number of bootstraps (b must = #filenames in i)\n"
					"       -c 			   		: number of K cluster in list k\n"
					"       -r             		: number of repetitions\n"
					"       -e filename    		: filename for storing cluster idx output after clusterMethodKmeans\n"
					"       -f filename    		: filename for storing cluster idx output after clusterMethodGmm\n"
					"       -g filename         : filename for storing cluster idx output after clusterMethodAgglo\n";

	fprintf(stderr, help, argv0);
}


double** readBootstrapData2D(char *bootstrapDataFiles, int *bootstrap_row, int *bootstrap_col) {

	int irow, icol;
	char *token;
	int i, j;
	size_t len;
	ssize_t numBytesRead;
	double **d_bootstrap_data2d;

	/* input file is in ASCII format -------------------------------*/
	FILE *infile;
	char *line, *ret;
	int lineLen;


						file_read_dim(bootstrapDataFiles,
								bootstrap_row, bootstrap_col);

						d_bootstrap_data2d = (double**) malloc(
								(*bootstrap_row) * sizeof(double*));
						assert(d_bootstrap_data2d != NULL);
						d_bootstrap_data2d[0] = (double*) malloc(
								(*bootstrap_row) * (*bootstrap_col) * sizeof(double));
						assert(d_bootstrap_data2d[0] != NULL);
						for (irow = 1; irow < (*bootstrap_row); irow++)
							d_bootstrap_data2d[irow] = d_bootstrap_data2d[irow
									- 1] + (*bootstrap_col);

						// reading file here itself as its giving so much problem in cluster
						//					printf("\nbootstrap_data1d malloc done..\n");

						/* first find the number of objects */
						lineLen = MAX_CHAR_PER_LINE;
						line = (char*) malloc(lineLen);

						if ((infile = fopen(strtok(bootstrapDataFiles, " \t\n"), "r")) == NULL) {

							d_printf("Error: no such file %s\n",
									bootstrapDataFiles);

							return NULL;
						}

						i = 0;
						while ((numBytesRead = getline(&line, &len, infile))
								!= -1) {
							token = strtok(line, " ");

							for (j = 0; j < (*bootstrap_col) && (token != NULL);
									j++) {

								d_bootstrap_data2d[i][j] = atof(token);

//								printf("%2.3f ", d_bootstrap_data2d[i][j]);

								token = strtok(NULL, " ");

							}
//							printf("\n");
							i++;
						}

						fclose(infile);
						free(line);

						return d_bootstrap_data2d;
}


double* readBootstrapData1D(char *bootstrapDataFiles, int *bootstrap_row, int *bootstrap_col) {

	int irow, icol;
	char *token;
	int i, j;
	size_t len;
	ssize_t numBytesRead;
	double *d_bootstrap_data1d;

	/* input file is in ASCII format -------------------------------*/
	FILE *infile;
	char *line, *ret;
	int lineLen;


						file_read_dim(bootstrapDataFiles,
								bootstrap_row, bootstrap_col);

						d_bootstrap_data1d = (double*) malloc((*bootstrap_row) * (*bootstrap_col) * sizeof(double));

						// reading file here itself as its giving so much problem in cluster
						//					printf("\nbootstrap_data1d malloc done..\n");

						/* first find the number of objects */
						lineLen = MAX_CHAR_PER_LINE;
						line = (char*) malloc(lineLen);


						if ((infile = fopen(strtok(bootstrapDataFiles, " \t\n"), "r")) == NULL) {

							d_printf("Error: no such file %s\n",
									bootstrapDataFiles);

							return NULL;
						}

						i = 0;
						while ((numBytesRead = getline(&line, &len, infile))
								!= -1) {
							token = strtok(line, " ");

							for (j = 0; j < (*bootstrap_col) && (token != NULL);
									j++) {

								d_bootstrap_data1d[j+i*(*bootstrap_col)] = atof(token);

//								printf("%2.3f ", d_bootstrap_data2d[i][j]);

								token = strtok(NULL, " ");

							}
//							printf("\n");
							i++;
						}

						fclose(infile);
						free(line);

						return d_bootstrap_data1d;
}



static void  malloc_init(void)
{
    mallopt(M_MMAP_MAX, 0);
}


int **alloc2d(int n, int m) {
	int i;
    int *data = (int*) malloc(n*m*sizeof(int));
    int **array = (int**) malloc(n*sizeof(int *));
    for (i=0; i<n; i++)
        array[i] = &(data[i*m]);
    return array;
}

int main(int argc, char** argv) {

	int opt;
	int master;
	char buf[READ_BUFF_LNGTH];
	/*extern char *optarg;*/
	extern int optind;
	int is_print_usage;
	_debug = 1;
	char *modelListFile;
	char *bootIdxFilesFile;
	char *modelList;
	char **bootstrapDataFiles;
	int bootstrapFileIdx;
	char *final_output_file;

	int nReps;
	int nK;

	FILE *ifp, *ofp;

	int partitionCount = 0;
	char *distmetric;
	char *centerfun;
	char *linkcode;
	char c_distmetric, c_centerfun;


	int count, irow, icol;

	int rank,spec_rank, nproc, mpi_namelen;
	char mpi_name[MPI_MAX_PROCESSOR_NAME];

	char *line;
	size_t len;
	ssize_t numBytesRead;
	int totalNumModels = 0;

	MPI_File fh;
	MPI_Status status;

	int *junk;
	int pool_size, number_of_blocks = 0;

	MPI_Offset my_offset, my_current_offset;

	MPI_Offset file_offset=0;

	char *const fmt="%4d ";
	char *const endfmt="%4d \n";
	const int charspernum=4;
	// modellist = 5, bootstrapFileIdx = 1, nK = 1, nReps = 1, timetaken = 3 (%4.4f)
	const int out_file_ncols = 11;
	MPI_Datatype num_as_string;
	MPI_Datatype localarray;

	char tempDataInputFile[TMP_FILE_NM_SIZE],
   		 tempDistOutputFile[TMP_FILE_NM_SIZE],
		 tempDist2SimFile[TMP_FILE_NM_SIZE],
		 tempEvdValsFile[TMP_FILE_NM_SIZE],
		 tempEvdVecFile[TMP_FILE_NM_SIZE],
		 tempClustIdxFile[TMP_FILE_NM_SIZE];

	int globalsizes[2];// = {totalNumModels, 5};
	int localsizes [2];// = {1, 5};
	int starts[2];//      = {0, 0};
	int order;//          = MPI_ORDER_C;
	int startrow, endrow, locnrows;

	int **idx;
	char *data_as_txt;


	int idx_ncols;

	MPI_Group orig_group, new_group;
	MPI_Comm new_comm;
	int *ranks_non_spectral, *ranks_spectral;

	double **d_bootstrap_data2d;
	double **temp_spec_bootstrap_data2d;

	int bootstrap_row, bootstrap_col, numBootstrap;

	// Variables to keep timings

	double model_start, model_end, global_start, global_end, local_start, local_end, cumulativeTimeTaken;
	double local_time;

	malloc_init();

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Get_processor_name(mpi_name, &mpi_namelen);

	/* some default values */
	_debug = 0;
	master = 0;
	is_print_usage = 0;
	nReps = 0;

	while ((opt = getopt(argc, argv, "m:o:n:f:b:s:r:c:")) != EOF)
		switch (opt) {
		case 'm':
			modelListFile = optarg;
			break;
		case 'n':
			totalNumModels = atoi(optarg);
			break;
		case 'o':
			final_output_file = optarg;
			break;
		case 'c':
			bootstrap_col  = atoi(optarg);
			break;
		case 'b':
			bootIdxFilesFile  = optarg;
			break;
		case 's':
			numBootstrap  = atoi(optarg);
			break;
		case 'r':
			bootstrap_row  = atoi(optarg);
			idx_ncols = bootstrap_row;
			break;

		default:
			is_print_usage = 1;
			break;
		}


	final_output_file = (char*) malloc(MAX_CHAR_PER_LINE);

	bootstrapDataFiles = (char**) malloc(numBootstrap * sizeof(char*));
	for(irow=0; irow<numBootstrap; irow++){
		bootstrapDataFiles[irow] = (char*) malloc(TMP_FILE_NM_SIZE * sizeof(char));
	}

	d_bootstrap_data2d = (double**) malloc(numBootstrap * sizeof(double*));
		assert(d_bootstrap_data2d != NULL);
		d_bootstrap_data2d[0] = (double*) malloc(bootstrap_row * bootstrap_col * numBootstrap * sizeof(double));
		assert(d_bootstrap_data2d[0] != NULL);
		for (irow = 1; irow < numBootstrap; irow++)
			d_bootstrap_data2d[irow] = d_bootstrap_data2d[irow - 1] + bootstrap_row*bootstrap_col;



	ifp = fopen(bootIdxFilesFile, "r");
	count = 0;

	while (fscanf(ifp, "%s", bootstrapDataFiles[count]) != EOF && count < numBootstrap) {
		d_printf("%s %d\n", bootstrapDataFiles[count], count);
		d_bootstrap_data2d[count] = readBootstrapData1D(bootstrapDataFiles[count],&bootstrap_row, &bootstrap_col);
		count++;
	}
	fclose(ifp);

//	for(count=0; count<numBootstrap; count++){
//		d_printf("%s %d\n", bootstrapDataFiles[count], count);
//		for(irow=0; irow<bootstrap_row; irow++){
//			for(icol=0;icol<bootstrap_col; icol++){
//				printf("%2.2f ", d_bootstrap_data2d[count][icol+irow*bootstrap_col]);
//			}
//			printf("\n");
//		}
//		printf("-----------------------------------------------------------------------------------\n");
//	}



//	line = (char*) malloc(MAX_CHAR_PER_LINE);
//
//	while ((numBytesRead = getline(&line, &len, ifp)) != -1) {
//
//		d_printf("%s\n", line);
//
//	}


    locnrows = totalNumModels/nproc;
    startrow = rank * locnrows;
    endrow = startrow + locnrows - 1;
    if (rank == nproc-1) {
        endrow = totalNumModels - 1;
        locnrows = endrow - startrow + 1;
    }

    ifp = fopen(modelListFile, "r");

	if(ifp==NULL){
		fprintf(stderr, "Error: no such file (%s)\n", modelListFile);
		return 0;

	}else{

	idx = alloc2d(totalNumModels, idx_ncols);


	/* convert our data into txt */
	data_as_txt = (char*) malloc(locnrows * (idx_ncols + out_file_ncols + 1) *charspernum*sizeof(char));
//	memset( data_as_txt, ' ', locnrows * (idx_ncols + out_file_ncols + 1) *charspernum * sizeof(char) );


	partitionCount = 0;
	count = 0;
	irow = 0;
	modelList = (char*) malloc(MAX_CHAR_PER_LINE);

		temp_spec_bootstrap_data2d = (double**) malloc(
				bootstrap_row * sizeof(double*));
		assert(temp_spec_bootstrap_data2d != NULL);
		temp_spec_bootstrap_data2d[0] = (double*) malloc(
				bootstrap_row * bootstrap_col * sizeof(double));
		assert(temp_spec_bootstrap_data2d[0] != NULL);
		for (irow = 1; irow < bootstrap_row; irow++)
			temp_spec_bootstrap_data2d[irow] = temp_spec_bootstrap_data2d[irow
					- 1] + bootstrap_col;

	local_start = MPI_Wtime();
	if(rank==0){
		global_start =  MPI_Wtime();
	}

	/* Extract the original group handle */
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	while (fscanf(ifp, "%s %d %d %d", modelList, &nK, &bootstrapFileIdx, &nReps) != EOF) {

		if (partitionCount++ % nproc != rank) {
		        continue;
		  }

		model_start = MPI_Wtime();

		distmetric = (char*) malloc(TMP_FILE_NM_SIZE * sizeof(char));
		centerfun = (char*) malloc(TMP_FILE_NM_SIZE * sizeof(char));
		linkcode = malloc(TMP_FILE_NM_SIZE * sizeof(char));


		// DO SPECTRAL
		if(getModelType(SPECTRAL_SHRT, modelList)){


			getDistMtrcnCntrFunBySpectral(SPECTRAL_SHRT, modelList, distmetric, centerfun);

			ranks_spectral = (int*) malloc(1 * sizeof(int));
			ranks_spectral[0] = rank;

			MPI_Group_incl(orig_group, 1, ranks_spectral, &new_group);
			/* Create new new communicator and then perform collective communications */
			MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);


			for(irow=0; irow< bootstrap_row; irow++){
				for(icol=0; icol<bootstrap_col; icol++){
					temp_spec_bootstrap_data2d[irow][icol] = d_bootstrap_data2d[bootstrapFileIdx][icol+irow*bootstrap_col];
				}
			}


			MPI_Group_rank (new_group, &spec_rank);


			spectral_adapter2D(new_comm, nK, bootstrap_row, bootstrap_col, temp_spec_bootstrap_data2d, centerfun, distmetric, idx[irow]);


			MPI_Comm_free(&new_comm);

			MPI_Group_free(&new_group);


		}
		// DO KMEANS
		else if(getModelType(KMEANS_SHRT, modelList))
		{
			getDistMtrcnCntrFunByKmeans(KMEANS_SHRT, modelList, distmetric, centerfun);

			for(irow=0; irow< bootstrap_row; irow++){
					for(icol=0; icol<bootstrap_col; icol++){
						temp_spec_bootstrap_data2d[irow][icol] = d_bootstrap_data2d[bootstrapFileIdx][icol+irow*bootstrap_col];
				}
			}

			kclust_adapter2D(nK, bootstrap_row,
											bootstrap_col, temp_spec_bootstrap_data2d, nReps,
											centerfun, distmetric, idx[irow]);


		}
		// DO AGGLOMERATIVE \ HIERARCHIAL
		else if(getModelType(AGGLO_SHRT, modelList))
		{
			getDistMtrcnCntrFunByAgglo(AGGLO_SHRT, modelList, distmetric, centerfun, linkcode);

			agglom_adapter1D(nK, bootstrap_row,
											bootstrap_col, d_bootstrap_data2d[bootstrapFileIdx],
											strtok(linkcode, " \t\n"), strtok(distmetric, " \t\n"), idx[irow]);

		}


		model_end = MPI_Wtime();

		// modellist width= 5*4
		sprintf(&data_as_txt[count*charspernum], "%s ", modelList);
//		printf(&data_as_txt[count*charspernum], "%s ", modelList);
		count+= 5;
		// bootstrapDataFiles width= 1*4
		sprintf(&data_as_txt[count*charspernum], "%4d ", bootstrapFileIdx);
//		printf(&data_as_txt[count*charspernum], "%4d ", bootstrapFileIdx);
		count+= 1;
		// nK width= 1*4
		sprintf(&data_as_txt[count*charspernum], "%4d ", nK);
//		printf(&data_as_txt[count*charspernum], "%4d ", nK);
		count++;
		// nReps width= 1*4
		sprintf(&data_as_txt[count*charspernum], "%4d ", nReps);
//		printf(&data_as_txt[count*charspernum], "%4d ", nReps);
		count++;

		// timetaken width= 2*4
		 sprintf(&data_as_txt[count*charspernum], "  %6.6f ", model_end - model_start);
		 count+=3;
		// increase width for next row
		for (icol=0; icol<idx_ncols; icol++) {
		    	sprintf(&data_as_txt[count*charspernum], fmt, idx[irow][icol]);
//		    	printf(&data_as_txt[count*charspernum], fmt, idx[irow][icol]);
		        count++;
		 }
		 sprintf(&data_as_txt[count*charspernum], "%s", "\n");
//		 printf(&data_as_txt[count*charspernum], "%s", "\n");
		 count++;

		free(distmetric);
		free(centerfun);

		irow++;


	}

	local_end = MPI_Wtime();
	local_time = local_end - local_start;
	MPI_Reduce(&local_time, &cumulativeTimeTaken, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	/* create a type describing our piece of the array */
	globalsizes[0] = totalNumModels;
	globalsizes[1] = idx_ncols + out_file_ncols;
	localsizes [0] = locnrows;
	localsizes [1] = idx_ncols + out_file_ncols;
	starts[0]      = startrow;
	starts[1]      = 0;
	order          = MPI_ORDER_C;

	/* each number is represented by charspernum chars */
    MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string);
    MPI_Type_commit(&num_as_string);


	MPI_File_open(MPI_COMM_WORLD, "output_idx.txt",
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_atomicity( fh, 1 ) ;


	MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order, num_as_string, &localarray);
	MPI_Type_commit(&localarray);


	MPI_File_set_view(fh, 0,  MPI_CHAR, localarray, "native", MPI_INFO_NULL);

	MPI_File_write_all(fh, data_as_txt, locnrows*idx_ncols, num_as_string, &status);

	MPI_File_close(&fh);

	close(&ifp);

	free(data_as_txt);
	free(idx);

	}



	MPI_Barrier(MPI_COMM_WORLD);



	d_printf("\n TIME TOOK RECORDED AT %d = %4.4f\n", rank, local_time);

	if(rank==0){
				global_end =  MPI_Wtime();
				printf("\n TIME TOOK RECORDED AT MASTER = %4.4f\n", global_end - global_start);
				printf("\n AVERAGE TIME TOOK RECORDED   = %4.4f\n", cumulativeTimeTaken/nproc);
	}

	MPI_Finalize();
	return 0;
}
