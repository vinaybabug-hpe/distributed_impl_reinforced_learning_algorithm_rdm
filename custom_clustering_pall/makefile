################################################################################
# Make file use to generate matlab executable custom ensemble clustering functions!
################################################################################

#-include ../makefile.init


RM := rm -rf
OUT = lib/myKmeans.a run_clust_one_mchine run_clust_p_mchine

LDIR=lib
ODIR = obj
SDIR = src
EXTERNAL=external
MEX_SDIR = mex_src
TESTS=tests
CPP_ODIR = obj/cpp
VLFEAT_DIR=external/vlfeat-0.9.20/bin/glnxa64



INC = -Iinc -Iexternal
CC = gcc
MEX_EXE=mex
MPICC= mpicc
MPICXX= mpicxx
CFLAGS=-fPIC -std=c99 -static -g -fno-inline -w

MPICC_CFLAGS=-g -fno-inline -fpermissive -Wwrite-strings -w

CFLAGS += -Wall -Wextra
CFLAGS += -Wno-unused-function -Wno-long-long -Wno-variadic-macros -Wsign-compare -Wunused-parameter -Wunused-variable

MEX_CFLAGS=-largeArrayDims -g 
MEX_LIBS=-Llib
DLIBS=-L$(LDIR) -Lexternal/parpack/ -Lexternal/f2c/
CPP_INC = -Iinc -Iexternal
CPP = g++
CPP_CFLAGS=-fPIC -std=c++0x -static -g -Wwrite-strings -O3 -Wall -Wno-sign-compare -fpermissive -fno-inline -w

# Select which features to disable for vlfeat
DISABLE_SSE2=yes
DISABLE_AVX=yes
DISABLE_THREADS=yes
DISABLE_OPENMP=yes

# All of the sources participating in the build are defined here

-include objects.mk

.PHONY: all clean pspectral

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: mkdir_o pspectral.lib mex2c_cluster_fit_pall build_seq_cluster_algo

# Tool invocations
kmeans_nwu: $(ODIR)/common/file_io.o $(ODIR)/nwu/seq_kmeans.o $(ODIR)/nwu/util.o $(ODIR)/main.o 
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc  -o "kmeans_nwu"  $(OBJS_NWU) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	
univ_of_tokyo_lib: $(ODIR)/common/file_io.o $(ODIR)/common/adapters.o $(ODIR)/univ_of_tokyo/cluster.o $(ODIR)/univ_of_tokyo/data.o  $(ODIR)/main.o  
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc  -o "univ_of_tokyo_lib"  $(OBJS_TOKYO) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	
mex2c_cluster_fit_pall: common_mex common 
	$(MEX_EXE) $(MEX_CFLAGS) CFLAGS="-fPIC -static" $(MEX_LIBS) $(INC) $(MEX_SDIR)/mex2c_cluster_fit_pall.c $(OBJS_COMMON_MEX) 

	
build_seq_cluster_algo: common vlfeat_static $(ODIR)/univ_of_tokyo/cluster.o $(ODIR)/univ_of_tokyo/data.o common_mex pspectral 
	$(MPICXX) $(MPICC_CFLAGS) $(INC) $(DLIBS) -o "seqClustAlgo" $(SDIR)/seqAlgoWrapper/seqClustAlgorithms.c $(OBJS_TOKYO_M) $(OBJS_COMMON_M) $(CPP_INC) $(LIBS) $(VLFEAT_LIB) $(OBJS_PSPECTRAL) $(LIBS) -lparpack -lf2c

build_pspec_cluster_algo: common  common_mex pspectral
	$(MPICXX) $(MPICC_CFLAGS) $(INC) $(DLIBS) -c $(ODIR)/pspectral/pSpectralOrchestrator.o $(OBJS_COMMON_MEX) $(OBJS_PSPECTRAL) $(CPP_INC) $(LIBS) -lparpack -lf2c  

#cluster_fit_threads:
#	$(CC) $(CFLAGS) -pthread $(INC) -o "cluster_fit_threads" $(SDIR)/cluster_fit_threads.c
	
	
vlfeat_static: 	$(ODIR)/vlfeat-0.9.20/vl/aib.o \
				$(ODIR)/vlfeat-0.9.20/vl/array.o \
				$(ODIR)/vlfeat-0.9.20/vl/covdet.o \
				$(ODIR)/vlfeat-0.9.20/vl/dsift.o \
				$(ODIR)/vlfeat-0.9.20/vl/fisher.o \
				$(ODIR)/vlfeat-0.9.20/vl/generic.o \
				$(ODIR)/vlfeat-0.9.20/vl/getopt_long.o \
				$(ODIR)/vlfeat-0.9.20/vl/gmm.o \
				$(ODIR)/vlfeat-0.9.20/vl/hikmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/homkermap.o \
				$(ODIR)/vlfeat-0.9.20/vl/host.o \
				$(ODIR)/vlfeat-0.9.20/vl/ikmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/imopv.o \
				$(ODIR)/vlfeat-0.9.20/vl/imopv_sse2.o \
				$(ODIR)/vlfeat-0.9.20/vl/kdtree.o \
				$(ODIR)/vlfeat-0.9.20/vl/kmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/lbp.o  \
				$(ODIR)/vlfeat-0.9.20/vl/liop.o \
				$(ODIR)/vlfeat-0.9.20/vl/mathop.o \
				$(ODIR)/vlfeat-0.9.20/vl/mathop_sse2.o \
				$(ODIR)/vlfeat-0.9.20/vl/mser.o \
				$(ODIR)/vlfeat-0.9.20/vl/pgm.o \
				$(ODIR)/vlfeat-0.9.20/vl/quickshift.o \
				$(ODIR)/vlfeat-0.9.20/vl/random.o \
				$(ODIR)/vlfeat-0.9.20/vl/rodrigues.o \
				$(ODIR)/vlfeat-0.9.20/vl/scalespace.o \
				$(ODIR)/vlfeat-0.9.20/vl/sift.o \
				$(ODIR)/vlfeat-0.9.20/vl/slic.o \
				$(ODIR)/vlfeat-0.9.20/vl/stringop.o \
				$(ODIR)/vlfeat-0.9.20/vl/svm.o \
				$(ODIR)/vlfeat-0.9.20/vl/svmdataset.o \
				$(ODIR)/vlfeat-0.9.20/vl/vlad.o \
				$(ODIR)/vlfeat/vlfeat_gmm.o
				ar  rcs $(LDIR)/libvl.a  \
				$(ODIR)/vlfeat-0.9.20/vl/aib.o \
				$(ODIR)/vlfeat-0.9.20/vl/array.o \
				$(ODIR)/vlfeat-0.9.20/vl/covdet.o \
				$(ODIR)/vlfeat-0.9.20/vl/dsift.o \
				$(ODIR)/vlfeat-0.9.20/vl/fisher.o \
				$(ODIR)/vlfeat-0.9.20/vl/generic.o \
				$(ODIR)/vlfeat-0.9.20/vl/getopt_long.o \
				$(ODIR)/vlfeat-0.9.20/vl/gmm.o \
				$(ODIR)/vlfeat-0.9.20/vl/hikmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/homkermap.o \
				$(ODIR)/vlfeat-0.9.20/vl/host.o \
				$(ODIR)/vlfeat-0.9.20/vl/ikmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/imopv.o \
				$(ODIR)/vlfeat-0.9.20/vl/imopv_sse2.o \
				$(ODIR)/vlfeat-0.9.20/vl/kdtree.o \
				$(ODIR)/vlfeat-0.9.20/vl/kmeans.o \
				$(ODIR)/vlfeat-0.9.20/vl/lbp.o  \
				$(ODIR)/vlfeat-0.9.20/vl/liop.o \
				$(ODIR)/vlfeat-0.9.20/vl/mathop.o \
				$(ODIR)/vlfeat-0.9.20/vl/mathop_sse2.o \
				$(ODIR)/vlfeat-0.9.20/vl/mser.o \
				$(ODIR)/vlfeat-0.9.20/vl/pgm.o \
				$(ODIR)/vlfeat-0.9.20/vl/quickshift.o \
				$(ODIR)/vlfeat-0.9.20/vl/random.o \
				$(ODIR)/vlfeat-0.9.20/vl/rodrigues.o \
				$(ODIR)/vlfeat-0.9.20/vl/scalespace.o \
				$(ODIR)/vlfeat-0.9.20/vl/sift.o \
				$(ODIR)/vlfeat-0.9.20/vl/slic.o \
				$(ODIR)/vlfeat-0.9.20/vl/stringop.o \
				$(ODIR)/vlfeat-0.9.20/vl/svm.o \
				$(ODIR)/vlfeat-0.9.20/vl/svmdataset.o \
				$(ODIR)/vlfeat-0.9.20/vl/vlad.o \
				$(ODIR)/vlfeat/vlfeat_gmm.o
				
				
# Tool invocations
c2cpp_test: mkdir_o cpp_wrapper_test
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	$(CC)  -o $(TESTS)/"c2cpp_wrapper_test" $(SDIR)/tests/C2CPPTest.c $(OBJS_C2CPP_TEST) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '				
				
common:	$(ODIR)/common/adapters.o $(ODIR)/common/utility_functions.o $(ODIR)/common/file_io.o
	
 
		 

common_mex: $(ODIR)/common/utility_functions.o $(ODIR)/common/file_io.o

cpp_wrapper_test: $(CPP_ODIR)/tests/c2cpp_cpptestfunc.o
	
pspectral.lib: 	parpacklib f2clib

pspectral: 	
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/pspectral/common.o $(SDIR)/pspectral/common.cc
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/pspectral/compute_distance.o $(SDIR)/pspectral/compute_distance.cc
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/pspectral/distance_to_similarity.o $(SDIR)/pspectral/distance_to_similarity.cc
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/pspectral/evd.o $(SDIR)/pspectral/evd.cc
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/pspectral/kmeans.o $(SDIR)/pspectral/kmeans.cc
	$(MPICXX) -c $(CPP_CFLAGS) $(CPP_INC) -o $(ODIR)/specwrapper/pSpectralWrapper.o $(SDIR)/specwrapper/pSpectralWrapper.cc
	
OBJS_PSPECTRAL= $(ODIR)/pspectral/common.o \
				$(ODIR)/pspectral/compute_distance.o \
				$(ODIR)/pspectral/distance_to_similarity.o \
				$(ODIR)/pspectral/evd.o \
				$(ODIR)/pspectral/kmeans.o \
				$(ODIR)/specwrapper/pSpectralWrapper.o
	
parpacklib:
	( cd external/parpack; $(MAKE))

f2clib:
	( cd external/f2c; $(MAKE))

#compile functions
_OBJS_COMMON = common/utility_functions.o common/file_io.o common/adapters.o
OBJS_COMMON_M = $(patsubst %,$(ODIR)/%,$(_OBJS_COMMON)) 

_OBJS_COMMON_MEX = common/utility_functions.o common/file_io.o 
OBJS_COMMON_MEX = $(patsubst %,$(ODIR)/%,$(_OBJS_COMMON_MEX)) 

_OBJS_NWU = common/file_io.o nwu/seq_kmeans.o nwu/util.o main.o
_OBJS_NWU_M = nwu/seq_kmeans.o nwu/util.o 
OBJS_NWU = $(patsubst %,$(ODIR)/%,$(_OBJS_NWU))
OBJS_NWU_M = $(patsubst %,$(ODIR)/%,$(_OBJS_NWU_M))

_OBJS_TOKYO = common/file_io.o common/adapters.o univ_of_tokyo/cluster.o univ_of_tokyo/data.o main.o
_OBJS_TOKYO_M = univ_of_tokyo/cluster.o univ_of_tokyo/data.o 
OBJS_TOKYO = $(patsubst %,$(ODIR)/%,$(_OBJS_TOKYO))
OBJS_TOKYO_M = $(patsubst %,$(ODIR)/%,$(_OBJS_TOKYO_M))

_OBJS_FILEIO = common/file_io.o
OBJS_FILEIO_M = $(patsubst %,$(ODIR)/%,$(_OBJS_FILEIO)) 

_OBJS_VLFEAT_M = vlfeat/vlfeat_gmm.o 
OBJS_VLFEAT_M = $(patsubst %,$(ODIR)/%,$(_OBJS_VLFEAT_M))

_OBJS_C2CPP_PSPEC = specwrapper/pSpectralWrapper.o specwrapper/pSpectralOrchestrator.o
OBJS_C2CPP_PSPEC = $(patsubst %,$(CPP_ODIR)/%,$(_OBJS_C2CPP_PSPEC))


#compile functions
_OBJS_C2CPP_TEST = tests/c2cpp_cpptestfunc.o
OBJS_C2CPP_TEST = $(patsubst %,$(CPP_ODIR)/%,$(_OBJS_C2CPP_TEST))


$(ODIR)/%.o: $(SDIR)/%.c mkdir_o
	$(CC) -c $(CFLAGS) $(INC) -o $@ $<  
	
$(CPP_ODIR)/%.o: $(SDIR)/%.cc  mkdir_o
	$(CPP) -c $(CPP_CFLAGS) $(CPP_INC) -o $@ $<  
	
mkdir_o:
	@mkdir -p $(ODIR)/common
	@mkdir -p $(ODIR)/nwu
	@mkdir -p $(ODIR)/univ_of_tokyo
	@mkdir -p $(ODIR)/vlfeat-0.9.20
	@mkdir -p $(ODIR)/vlfeat-0.9.20/vl
	@mkdir -p $(ODIR)/vlfeat
	@mkdir -p $(ODIR)/pspectral
	@mkdir -p $(CPP_ODIR)
	@mkdir -p $(CPP_ODIR)/lib_clustering
	@mkdir -p $(CPP_ODIR)/specwrapper
	@mkdir -p $(CPP_ODIR)/$(TESTS)
	@mkdir -p $(TESTS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*/*.o  $(CPP_ODIR)/*/*.o $(OUT) $(ODIR)/vlfeat-0.9.20/*/*.o $(ODIR)/*.o kmeans_nwu *.mexa64 seqClustAlgo pSpectralOrchestrator
	rm -f *.o external/parpack/*.o external/parpack/*.a external/f2c/*.o external/f2c/*.a
	
