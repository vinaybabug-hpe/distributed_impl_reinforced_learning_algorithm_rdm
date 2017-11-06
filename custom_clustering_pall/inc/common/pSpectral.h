/*
 ============================================================================
 Name        : pSpectral.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Oct 31, 2015
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


#ifndef PSPECTRAL_H_
#define PSPECTRAL_H_

/*
void pSpectralOrchestrator(double **data2d, int nrows, int ncols, char dist,
		int *clusterid, int nclusters, int numThreads, char* tempDataInputFile,
		char* tempDistOutputFile, char* tempDist2SimFile, char* tempEvdValsFile,
		char* tempEvdVecFile, char* tempClustIdxFile);
*/



void spectral_adapter2D(MPI_Comm new_comm, int nclusters, int nrows, int ncols, double** data2d,
		char* _method, char* _dist, int *clusterid/*, char* tempDataInputFile,
		char* tempDistOutputFile, char* tempDist2SimFile, char* tempEvdValsFile,
		char* tempEvdVecFile, char* tempClustIdxFile*/);

void internalOrchestrator(MPI_Comm new_comm, double **data2d, int nrows, int ncols, char dist,
		int *clusterid, int nclusters/*, char* tempDataInputFile,
		char* tempDistOutputFile, char* tempDist2SimFile, char* tempEvdValsFile,
		char* tempEvdVecFile, char* tempClustIdxFile*/);



#endif /* PSPECTRAL_H_ */
