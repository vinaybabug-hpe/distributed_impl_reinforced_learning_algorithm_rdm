/*
 ============================================================================
 Name        : wrapperFunc.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Oct 8, 2015
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

#ifndef INC_COMMON_WRAPPERFUNC_H_
#define INC_COMMON_WRAPPERFUNC_H_



int createTempFile(char *fileName);
int getModelType(char *clustMethod, char *methodList);
void getClusteringMethodsList(int length, char *methodList[], char **clustList);
int doKmeans(char **clustList);
int doKmedoids(char **clustList);
int doGMM(char **clustList);
int doSpectral(char **clustList);
int doAgglom(char **clustList);
int getMthdLstbyClust(char*clustMethod, char **methodList,
		int length, char **mthdLstbyClust);
int getClustMthdCnt(char*clustMethod, char **methodList,
		int length);
void getDistMtrcnCntrFunByKmeans(char*clustMthd, char *model,
		char *distmetric, char *centerfun);
void getDistMtrcnCntrFunByAgglo(const char*clustMthd, const char *model,
		char *distmetric, char *centerfun, char *linkcode);
void getDistMtrcnCntrFunBySpectral(const char*clustMthd, const char *model,
		char *distmetric, char *centerfun);

int kmeans_adapter(int nclusters, int nrow, int ncol, double* data1d, int npass,
		char method, char dist, int clusterid[]);
int kclust_adapter2D(int nclusters, int nrow, int ncol, double** data2d,
		int npass, char* _method, char* _dist, int clusterid[]);
void agglom_adapter1D(int nclusters, int nrows, int ncols, double* data1d,
		char* _method, char* _dist, int *clusterid);
void agglom_adapter2D(int nclusters, int nrows, int ncols, double** data2d,
		char* _method, char* _dist, int *clusterid);
void agglom_adapter1D(int nclusters, int nrows, int ncols, double* data1d,
		char* _method, char* _dist, int *clusterid);

void agglom_adapter(int nclusters, int nrows, int ncols, double* data1d,
		char method, char dist, int clusterid[]);
void kmedoids_adapter(int nclusters, int npass, int nrows, int ncols,
		double* data1d, char method, char dist, int clusterid[]);
void cluster_util_bootpartition2partition(unsigned short *boot, float *bootIdxs,
		int *idx, int length);



void vlfeat_gmm_adapter(int nclusters, int nrows,
		int ncols, double* data, int clusterid[]);

/* ========================================================================= */

float** file_read_ascii(
                  char *filename,      /* input file name */
                  int  *numObjs,       /* no. data objects (local) */
                  int  *numCoords)     /* no. coordinates */;


#endif /* INC_COMMON_WRAPPERFUNC_H_ */
