/*
 ============================================================================
 Name        : wrapper.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Aug 16, 2015
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

#ifndef WRAPPER_H_
#define WRAPPER_H_
#include <stdbool.h>

#define PRINT_DEBUG		0
#define TMP_NKLIST 		2
#define TMP_NREPS 		1

#define RUN_KMEANS 		0
#define RUN_AGGLO 		0
#define RUN_GMM		 	0
#define RUN_SPECTRAL 	1
#define CLUST_IN_ENSEMBLE 3

#define OUTPUT_FILE "data/output.txt"
#define KMEAN_TIMESTAMP_FILE "kmeans_run_times.txt"
#define AGGLO_TIMESTAMP_FILE "agglo_run_times.txt"
#define SPEC_TIMESTAMP_FILE "spec_run_times.txt"
#define GMM_TIMESTAMP_FILE "gmm_run_times.txt"
#define SEQ_CLUST_ALGO_EXECUTABLE "../seqClustAlgo"
#define PSPEC_CLUST_ALGO_EXECUTABLE "../pSpectralOrchestrator"
#define TESTCASE1_SH "testcase1.sh"


#define d_printf(fmt, args...) if(PRINT_DEBUG) { printf("[%s | %d | %s]: ",__FILE__,  __LINE__,__FUNCTION__); printf(fmt, ## args ); }

/**
 * % gmm, kmeansxxx, medoidxxx, spectralxxx and aggxxxyyy  where:
 * %  xxx denotes distance metric  (euc,seu,cit,cor,cos,mah,che,spe,ham,jac)
 * %      (for kmeans, can only use: (euc,cit,cor,cos,ham)
 * %   yyy denotes linkage metric     (avg,cen,com,med,sin,war,wei)
 */
#define KMEANS_SHRT 		"kme"
#define KMEDOIDS_SHRT 		"med"
#define GMM_SHRT 			"gmm"
#define SPECTRAL_SHRT 		"spe"
#define AGGLO_SHRT 			"agg"

#define KMEANS_LNG 			"kmeans"
#define KMEDOIDS_LNG 		"medoid"
#define GMM_LNG 			"gmm"
#define SPECTRAL_LNG 		"spectral"
#define AGGLO_LNG 			"agg"

#define DIST_MTRC_EUC		"euc" // euclidean
#define DIST_MTRC_SEU		"seu" // seuclidean
#define DIST_MTRC_CIT		"cit" // cityblock
#define DIST_MTRC_COR		"cor" // correlation
#define DIST_MTRC_COS		"cos" // cosine
#define DIST_MTRC_HAM		"ham" // hamming
#define DIST_MTRC_MAH		"mah" // mahalanobis
#define DIST_MTRC_JAC		"jac" // jaccard
#define DIST_MTRC_CHE		"che" // chebychev
#define DIST_MTRC_SPE		"spe" // spearman
#define DIST_MTRC_CODE_LN		3

#define CNTR_FUN_MEAN 		"mean"
#define CNTR_FUN_MEDIAN 	"median"

#define LNK_CODE_AVG	"avg" // average
#define LNK_CODE_CEN	"cen" // centroid
#define LNK_CODE_COM	"com" // complete
#define LNK_CODE_SIN	"sin" // single
#define LNK_CODE_MED	"med" // median
#define LNK_CODE_WAR	"war" // ward
#define LNK_CODE_WEI	"wei" // weighted
#define LNK_CODE_LN		3

/*#define TMP_FILE_FMT "/tmp/bootIdxTmpFile-XXXXXX"*/

#define TMP_FILE_FMT "tmp/bootIdxTmpFile-XXXXXX"

#define MAX_CHAR_PER_LINE 1024
#define TMP_FILE_NM_SIZE	64
#define TMP_FILE_MODEL_LIST "/tmp/ensemble_clust_mdl_lst"
#define CONFIG_FILE "config.cfg"
#define READ_BUFF_LNGTH		1000
#define DEBUG 0

#define NUM_MPI_THREADS     1
#define MPIRUN_CLUST_SING_MCHINE_EXE "run_clust_p_mchine"

/*
 * Field names for PST passed by RDM Ensemble Clustering
 */
#define NAME 							"name"
#define VERBOSE							"verboseQ"
#define OUTLIER_CUTOFF 					"outlierCutoff"
#define OUTLIER_METRIC 					"outlierMetric"
#define PREPROCESS_LIST 				"preprocessList"
#define ENSEMBLE						"ensemble"
#define EXTRACT 						"extract"
#define VALIDATE 						"validate"
#define COMPARE 						"compare"
#define KMEANS 							"kmeans"
#define MEDOID 							"medoid"
#define GMM 							"gmm"
#define SPECTRAL 						"spectral"
#define AGGLO 							"agg"

#define ENSEMBLE_KLIST 					"kList"
#define ENSEMBLE_MODELLIST 				"modelList"
#define ENSEMBLE_NBOOTSTRAPS 			"nBootstraps"
#define ENSEMBLE_NREPS					"nReps"
#define ENSEMBLE_INCLUDEREPSQ			"includeRepsQ"
#define ENSEMBLE_INCLUDECENTERSQ		"includeCentersQ"
#define ENSEMBLE_VERBOSEQ				"verboseQ"
#define MATLAB_BOOTIDX_FUNC				"getBootstrapData"


#define NaN								9999

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))
#define NUMBER_OF_FIELDS_PARTITIONS (sizeof(partitions_field_names)/sizeof(*partitions_field_names))
#define NUMBER_OF_FIELDS_CENTERS (sizeof(centers_field_names)/sizeof(*centers_field_names))



#endif /* WRAPPER_H_ */
