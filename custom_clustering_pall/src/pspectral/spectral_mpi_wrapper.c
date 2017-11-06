/*
 * spectral_wrapper.c
 *
 *  Created on: Aug 15, 2015
 *      Author: Vinay B Gavirangaswamy
 */

#include <stdlib.h> /* needed to define exit() */
#include <unistd.h> /* needed to define getpid() */
#include <stdio.h> /* needed for printf() */
#include <sys/types.h>
#include <string.h>
#include <sys/errno.h>
#include <fcntl.h>

#include "wrapper.h"



int createTempFile(char *fileName) {
	char nameBuff[32];
	int filedes = -1, errno_ = 0;
	// memset the buffers to 0
	memset(nameBuff, 0, sizeof(nameBuff));
	// Copy the relevant information in the buffers
	strncpy(nameBuff, TMP_FILE_FMT, 32);
	//string tempfile(nameBuff);

	// Create the temporary file, for input data
	filedes = mkstemp(nameBuff);
	strcpy(fileName, nameBuff);

	if (filedes < 1) {
		printf("\n Creation of temp file failed with error [%s]\n",
				strerror(errno_));
		return 1;
	} else {
		printf("\n Temporary file [%s] created\n", nameBuff);
	}
}

void printFile(char *tempClustIdxFile)
{
	char buf[1000];
	FILE *ifp;
    ifp = fopen(tempClustIdxFile, "r");
    while(fgets(buf, 1000, ifp) != NULL){
        // Each line is a training document.
        if(sizeof (buf) > 0 && // Skip empty lines
        		buf[0] != '\r' && // Skip empty lines
        		buf[0] != '\n' && // Skip empty lines
        		buf[0] != '#'){
            // Skip comment lines.
            //cout << line << endl;
            printf("%s", buf);
        }
    }

    fclose(ifp);
}

int main(int argc, char **argv) {

	int NUM_CLUSTERS = 10;
	char *inputFifoName="testdata/data.txt";
	char* tempClustIdxFile ="testdata/clusterIdxs.txt";;


	char buf[1000];
	int nbytes;
	pid_t childpid;
	char tempDataInputFile[32],  tempDistOutputFile[32], tempDist2SimFile[32], tempEvdValsFile[32], tempEvdVecFile[32];
	FILE *ifp, *ofp;


	int  errno_ = 0, start = START;

	createTempFile(tempDataInputFile);
	createTempFile(tempDistOutputFile);
	createTempFile(tempDist2SimFile);
	createTempFile(tempEvdValsFile);
	createTempFile(tempEvdVecFile);
	createTempFile(tempClustIdxFile);





	ifp = fopen(inputFifoName, "r");
	ofp = fopen(tempDataInputFile, "w");

	  int count = 0;
	  while (fgets(buf,1000, ifp)!=NULL) {  // Each line is a training document.
	    if (sizeof(buf) > 0 &&      // Skip empty lines.
	    	buf[0] != '\r' &&      // Skip empty lines.
	    	buf[0] != '\n' &&      // Skip empty lines.
	    	buf[0] != '#') {       // Skip comment lines.
	    	//cout << line << endl;
	    	count++;
	    	printf("%s", buf);
	    	fprintf(ofp,"%s", buf);
	    }
	  }

	  fclose(ifp);
	  fclose(ofp);
	  //unlink(tempInputFile);


	/*Step 1: Compute Similarity*/
	/*Spawn a child to run the program.*/
	childpid = fork();
	if (childpid == 0) { /* child process */
		//char *av[] = { mpiReadFd, mpiWriteFd, (char *) 0 }; /* each element represents a command line argument */

		//execv("Matlab_MPI", av);

		char **args = (char **) calloc(10, sizeof(char *));
		args[0] = "mpiexec";
		args[1] = "-n";
		args[2] = "8";
		args[3] = "compute_distance";
		args[4] = "--t_nearest_neighbor";
		args[5] = "8"; // input argument
		args[6] = "--input";
		args[7] = tempDataInputFile;
		args[8] = "--output";
		args[9] = tempDistOutputFile;

		//args[6] = (char *) 0;

		execvp("mpiexec", args);
		perror("execvp"); /* if we get here, execvp failed */

	} else { /* pid!=0; parent process */
		/*printf("\nParent continues..\n");*/
		waitpid(childpid, 0, 0); /* wait for child to exit */
	}
	printf("\nStep 1: Compute Distance \t[DONE]\n");

	/*Step 2: Compute distance_to_similarity*/
	/*Spawn a child to run the program.*/
	childpid = fork();
	if (childpid == 0) { /* child process */
		//char *av[] = { mpiReadFd, mpiWriteFd, (char *) 0 }; /* each element represents a command line argument */

		//execv("Matlab_MPI", av);

		char **args = (char **) calloc(8, sizeof(char *));
		args[0] = "mpiexec";
		args[1] = "-n";
		args[2] = "8";
		args[3] = "distance_to_similarity";
		args[4] = "--input";
		args[5] = tempDistOutputFile;
		args[6] = "--output";
		args[7] = tempDist2SimFile;

		//args[6] = (char *) 0;

		execvp("mpiexec", args);
		perror("execvp"); /* if we get here, execvp failed */

	} else { /* pid!=0; parent process */
		/*printf("\nParent continues..\n");*/
		waitpid(childpid, 0, 0); /* wait for child to exit */
	}
	printf("\nStep 2: Compute distance_to_similarity \t[DONE]\n");
	printFile(tempDist2SimFile);

	/*Step 3: Compute evd*/
	/*Spawn a child to run the program.*/
	childpid = fork();
	if (childpid == 0) { /* child process */
		//char *av[] = { mpiReadFd, mpiWriteFd, (char *) 0 }; /* each element represents a command line argument */

		//execv("Matlab_MPI", av);

		char **args = (char **) calloc(14, sizeof(char *));
		args[0] = "mpiexec";
		args[1] = "-n";
		args[2] = "16";
		args[3] = "evd";
		args[4] = "--eigenvalue";
		args[5] = "100";
		args[6] = "--eigenspace";
		args[7] = "300";
		args[8] = "--input";
		args[9] = tempDist2SimFile;
		args[10] = "--eigenvalues_output";
		args[11] = tempEvdValsFile;
		args[12] = "--eigenvectors_output";
		args[13] = tempEvdVecFile;

		//args[6] = (char *) 0;

		execvp("mpiexec", args);
		perror("execvp"); /* if we get here, execvp failed */

	} else { /* pid!=0; parent process */
		/*printf("\nParent continues..\n");*/
		waitpid(childpid, 0, 0); /* wait for child to exit */
	}
	printf("\nStep 3: Compute evd \t[DONE]\n");

	/*Step 4: Compute kmeans on evd for idx*/
	/*Spawn a child to run the program.*/
	childpid = fork();
	if (childpid == 0) { /* child process */
		//char *av[] = { mpiReadFd, mpiWriteFd, (char *) 0 }; /* each element represents a command line argument */

		//execv("Matlab_MPI", av);

		sprintf(buf, "%d", NUM_CLUSTERS);

		char **args = (char **) calloc(10, sizeof(char *));
		args[0] = "mpiexec";
		args[1] = "-n";
		args[2] = "8";
		args[3] = "kmeans";
		args[4] = "--num_clusters";
		args[5] = buf;
		args[6] = "--input";
		args[7] = tempEvdVecFile;
		args[8] = "--output";
		args[9] = tempClustIdxFile;


		//args[6] = (char *) 0;

		execvp("mpiexec", args);
		perror("execvp"); /* if we get here, execvp failed */

	} else { /* pid!=0; parent process */
		/*printf("\nParent continues..\n");*/
		waitpid(childpid, 0, 0); /* wait for child to exit */
	}
	printf("\nStep 4: Compute kmeans on evd for idx \t[DONE]\n");
    printFile(tempClustIdxFile);

	unlink(tempDataInputFile);
	unlink(tempDistOutputFile);
	unlink(tempDist2SimFile);
	unlink(tempEvdValsFile);
	unlink(tempEvdVecFile);
	//unlink(tempClustIdxFile);
	printf("EXECUTED!\n");

	exit(1);
}
