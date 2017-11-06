#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <cassert>
#include <Eigen/Core>
#include "SpectralClustering.h"
#include "spectral_wrapper.h"

extern "C" { // another way
float** file_read_ascii(char *filename, /* input file name */
int *numObjs, /* no. data objects (local) */
int *numCoords); /* no. coordinates */

double** getDistMetric(int nrows, int ncols, double** data, char dist);

};




int main() {

	char fileName[32] = "color1002.txt";
	int data_row, data_col;
	float **data;
	int i, j;
	char dist = 'e';
	double **distmatrix;
	double sigma = 0.04;
	int *idx;


	data = file_read_ascii(fileName, &data_row, &data_col);
	//centralize and scale the data

	double **data_double;

	idx = (int*)malloc(data_row * sizeof(int));

	data_double = create2dMat(data_row, data_col);

	for (i = 0; i < data_row; i++) {
		for (j = 0; j < data_col; j++) {
			data_double[i][j] = data[i][j];
		}
	}

//	normalize_data(data_double, data_row, data_col);
	distmatrix =getDistMetric(data_row, data_col, data_double, dist);

//	distMatrix = create2dMat(data_row, data_row);
//
//
//
//	if(temp_distmatrix != NULL)
//	for (i = 0; i < data_row; i++) {
//		for (j = 0; j < i; j++) {
//			//std::cout <<" distmatrix[i][j]= "<< temp_distmatrix[i][j] << " ";
//			distMatrix[i][j]=temp_distmatrix[i][j];
//			distMatrix[j][i]=temp_distmatrix[i][j];
//		}
//		distMatrix[i][i]=0;
//	}



	std::cout << "Starting Spectal Clustering!" << std::endl;

    // generate similarity matrix
    unsigned int size = data_row;
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size,size);

    for (unsigned int i=0; i < size; i++) {
        for (unsigned int j=0; j < i; j++) {
            // generate similarity

            double similarity = exp(-distmatrix[i][j]*distmatrix[i][j] / (2*sigma*sigma));
            m(i,j) = similarity;
            m(j,i) = similarity;
           // std::cout<<m(i,j)<<" ";
        }
//        std::cout<<std::endl;
        m(i,i) = 0;
    }

    //std::cout << "generate similarity matrix [DONE!]" << std::endl;

    // the number of eigenvectors to consider. This should be near (but greater) than the number of clusters you expect. Fewer dimensions will speed up the clustering
    int numDims = size;

    // do eigenvalue decomposition
    SpectralClustering* c = new SpectralClustering(m, numDims);

    // whether to use auto-tuning spectral clustering or kmeans spectral clustering
    bool autotune = false;

    std::vector<std::vector<int> > clusters;
    if (autotune) {
        // auto-tuning clustering
        clusters = c->clusterRotate();
    } else {
        // how many clusters you want
        int numClusters = 5;
        clusters = c->clusterKmeans(numClusters);
    }


    for (unsigned int i=0; i < clusters.size(); i++) {
    	std::cout << "Cluster " << i << ": " << "Item ";
        std::copy(clusters[i].begin(), clusters[i].end(), std::ostream_iterator<int>(std::cout, ", "));
        std::cout << std::endl;
    }

    // output clustered items
    // items are ordered according to distance from cluster centre
    for (unsigned int i=0; i < clusters.size(); i++) {
    	int temp;
    	temp = clusters[i].size();
    	for(j=0; j<temp; j++){

    		idx[clusters[i].at(j)] = i;

    	}

    }

    for(i=0; i< data_row; i++){
    	std::cout <<"["<< i<<"]"<<" = " <<idx[i]<<std::endl;
    }

    std::cout<<std::endl;
}
