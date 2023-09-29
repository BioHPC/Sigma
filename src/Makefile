############
# Makefile #
############


#########################################################
# Please change values based on your system environment #
#########################################################

# C++ Compiler
CC = g++

# Open MPI Compiler
MPICC = mpiCC

# C++ Compile Flag
CFLAGS = -c -Wall

# GCC OpenMP Flag
OMPFLAG = -fopenmp

# Static build
STFLAG = -static

# Include directories 
INCL = -I./Ipopt/include

# Library
#LDFLAGS = -L./Ipopt/lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../x86_64-suse-linux/lib -L/usr/lib64/gcc/x86_64-suse-linux/4.7/../../.. -lipopt -llapack -ldl -lcoinmumps -lfblas -lpthread -lgfortran -lm -lquadmath
#LDFLAGS = -L./Ipopt/lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../.. -lipopt -llapack -ldl -lcoinmumps -lblas -lgfortran -lm -lquadmath
LDFLAGS = -L./Ipopt/lib64 -L/lib64 -lipopt -llapack -ldl -lcoinmumps -lblas -lgfortran -lm -lquadmath


#########################################################
# Please do NOT change below Makefile variables         #
#########################################################

# Common source files
SIGMA_COMMON_SRC = utils.cpp \
	parse_config.cpp \
	sigma_config.cpp \
	sigma_core.cpp \
	sigma_stat.cpp

# SIGMA index genomes
SIGMA_INDEX_GENOMES_SRC = sigma_index_genomes_main.cpp 
SIGMA_INDEX_GENOMES_MPI_SRC = sigma_index_genomes_mpi_main.cpp 

# SIGMA align reads
SIGMA_ALIGN_READS_SRC = sigma_align_reads_main.cpp 
SIGMA_ALIGN_READS_MPI_SRC = sigma_align_reads_mpi_main.cpp 

# SIGMA build model
SIGMA_BUILD_MODEL_SRC = sigma_build_model_main.cpp

# SIGMA solve model
SIGMA_SOLVE_MODEL_SRC = sigma_solve_model_main.cpp \
	sigma_ipopt_run.cpp \
	sigma_ipopt_nlp.cpp

# SIGMA main
SIGMA_SRC = sigma_main.cpp \
	sigma_ipopt_run.cpp \
	sigma_ipopt_nlp.cpp

# SIGMA stat bootstrap
SIGMA_BOOTSTRAP_MPI_SRC = sigma_bootstrap_mpi_main.cpp 

# SIGMA stat bootstrap
SIGMA_JACKKNIFE_MPI_SRC = sigma_jackknife_mpi_main.cpp 

# Object files are from source files
SIGMA_COMMON_OBJS = ${SIGMA_COMMON_SRC:.cpp=.o}
SIGMA_INDEX_GENOMES_OBJS = ${SIGMA_INDEX_GENOMES_SRC:.cpp=.o}
SIGMA_INDEX_GENOMES_MPI_OBJS = ${SIGMA_INDEX_GENOMES_MPI_SRC:.cpp=.o}
SIGMA_ALIGN_READS_OBJS = ${SIGMA_ALIGN_READS_SRC:.cpp=.o}
SIGMA_ALIGN_READS_MPI_OBJS = ${SIGMA_ALIGN_READS_MPI_SRC:.cpp=.o}
SIGMA_BUILD_MODEL_OBJS = ${SIGMA_BUILD_MODEL_SRC:.cpp=.o}
SIGMA_SOLVE_MODEL_OBJS = ${SIGMA_SOLVE_MODEL_SRC:.cpp=.o}
SIGMA_OBJS = ${SIGMA_SRC:.cpp=.o}
SIGMA_BOOTSTRAP_MPI_OBJS = ${SIGMA_BOOTSTRAP_MPI_SRC:.cpp=.o}
SIGMA_JACKKNIFE_MPI_OBJS = ${SIGMA_JACKKNIFE_MPI_SRC:.cpp=.o}

# This should be the name of your executable
SIGMA_INDEX_GENOMES_EXE = sigma-index-genomes
SIGMA_INDEX_GENOMES_MPI_EXE = sigma-index-genomes-mpi
SIGMA_ALIGN_READS_EXE = sigma-align-reads
SIGMA_ALIGN_READS_MPI_EXE = sigma-align-reads-mpi
SIGMA_BUILD_MODEL_EXE = sigma-build-model
SIGMA_SOLVE_MODEL_EXE = sigma-solve-model
SIGMA_EXE = sigma
SIGMA_BOOTSTRAP_MPI_EXE = sigma-bootstrap-mpi
SIGMA_JACKKNIFE_MPI_EXE = sigma-jackknife-mpi

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

.SUFFIXES: .cpp .c .o .obj

.PHONY: all mpi

# Build arg

all: sigma_index_genomes sigma_align_reads sigma_build_model sigma_solve_model sigma

mpi: sigma_index_genomes_mpi sigma_align_reads_mpi sigma_bootstrap_mpi sigma_jackknife_mpi

sigma_index_genomes: $(SIGMA_INDEX_GENOMES_EXE)

sigma_index_genomes_mpi: $(SIGMA_INDEX_GENOMES_MPI_EXE)

sigma_align_reads: $(SIGMA_ALIGN_READS_EXE)

sigma_align_reads_mpi: $(SIGMA_ALIGN_READS_MPI_EXE)

sigma_build_model: $(SIGMA_BUILD_MODEL_EXE)

sigma_solve_model: $(SIGMA_SOLVE_MODEL_EXE)

sigma: $(SIGMA_EXE)

sigma_bootstrap_mpi: $(SIGMA_BOOTSTRAP_MPI_EXE)

sigma_jackknife_mpi: $(SIGMA_JACKKNIFE_MPI_EXE)


# Build

$(SIGMA_INDEX_GENOMES_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_INDEX_GENOMES_OBJS)
	$(CC) $^ -o $@

$(SIGMA_INDEX_GENOMES_MPI_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_INDEX_GENOMES_MPI_OBJS)
	$(MPICC) $^ -o $@

$(SIGMA_ALIGN_READS_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_ALIGN_READS_OBJS)
	$(CC) $^ -o $@

$(SIGMA_ALIGN_READS_MPI_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_ALIGN_READS_MPI_OBJS)
	$(MPICC) $^ -o $@

$(SIGMA_BUILD_MODEL_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_BUILD_MODEL_OBJS)
	$(CC) $^ -o $@

$(SIGMA_SOLVE_MODEL_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_SOLVE_MODEL_OBJS)
	$(CC) $(OMPFLAG) $(STFLAG) $^ $(LDFLAGS) -o $@

$(SIGMA_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_OBJS)
	$(CC) $(OMPFLAG) $(STFLAG) $^ $(LDFLAGS) -o $@

$(SIGMA_BOOTSTRAP_MPI_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_BOOTSTRAP_MPI_OBJS)
	$(MPICC) $^ -o $@

$(SIGMA_JACKKNIFE_MPI_EXE): $(SIGMA_COMMON_OBJS) $(SIGMA_JACKKNIFE_MPI_OBJS)
	$(MPICC) $^ -o $@


# Clean

clean:
	rm -f $(SIGMA_INDEX_GENOMES_EXE) \
	$(SIGMA_INDEX_GENOMES_MPI_EXE) \
	$(SIGMA_ALIGN_READS_EXE) \
	$(SIGMA_ALIGN_READS_MPI_EXE) \
	$(SIGMA_BUILD_MODEL_EXE) \
	$(SIGMA_SOLVE_MODEL_EXE) \
	$(SIGMA_EXE) \
	$(SIGMA_BOOTSTRAP_MPI_EXE) \
	$(SIGMA_JACKKNIFE_MPI_EXE) \
	$(SIGMA_COMMON_OBJS) \
	$(SIGMA_INDEX_GENOMES_OBJS)	\
	$(SIGMA_INDEX_GENOMES_MPI_OBJS)	\
	$(SIGMA_ALIGN_READS_OBJS) \
	$(SIGMA_ALIGN_READS_MPI_OBJS) \
	$(SIGMA_BUILD_MODEL_OBJS) \
	$(SIGMA_SOLVE_MODEL_OBJS) \
	$(SIGMA_OBJS) \
	$(SIGMA_BOOTSTRAP_MPI_OBJS) \
	$(SIGMA_JACKKNIFE_MPI_OBJS)
	
####
# To handle CC and MPICC, let compile individually
#.cpp.o:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) -o $@ $<
#
#.cpp.obj:
#    $(CC) $(CFLAGS) $(EXTRA_FLAGS) $(INCL) $^ -o $@ `$(CYGPATH_W) '$<'`

utils.o: utils.cpp
	$(CC) $(CFLAGS) utils.cpp

parse_config.o: parse_config.cpp
	$(CC) $(CFLAGS) parse_config.cpp

sigma_config.o: sigma_config.cpp
	$(CC) $(CFLAGS) sigma_config.cpp

sigma_core.o: sigma_core.cpp
	$(CC) $(CFLAGS) sigma_core.cpp

sigma_stat.o: sigma_stat.cpp
	$(CC) $(CFLAGS) sigma_stat.cpp

$(SIGMA_INDEX_GENOMES_OBJS): $(SIGMA_INDEX_GENOMES_SRC)
	$(CC) $(CFLAGS) $(SIGMA_INDEX_GENOMES_SRC)

$(SIGMA_INDEX_GENOMES_MPI_OBJS): $(SIGMA_INDEX_GENOMES_MPI_SRC)
	$(MPICC) $(CFLAGS) $(SIGMA_INDEX_GENOMES_MPI_SRC)

$(SIGMA_ALIGN_READS_OBJS): $(SIGMA_ALIGN_READS_SRC)
	$(CC) $(CFLAGS) $(SIGMA_ALIGN_READS_SRC)

$(SIGMA_ALIGN_READS_MPI_OBJS): $(SIGMA_ALIGN_READS_MPI_SRC)
	$(MPICC) $(CFLAGS) $(SIGMA_ALIGN_READS_MPI_SRC)

$(SIGMA_BUILD_MODEL_OBJS): $(SIGMA_BUILD_MODEL_SRC)
	$(CC) $(CFLAGS) $(SIGMA_BUILD_MODEL_SRC)

sigma_solve_model_main.o: sigma_solve_model_main.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_solve_model_main.cpp

sigma_ipopt_run.o: sigma_ipopt_run.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_ipopt_run.cpp

sigma_ipopt_nlp.o: sigma_ipopt_nlp.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_ipopt_nlp.cpp

sigma_main.o: sigma_main.cpp
	$(CC) $(CFLAGS) $(OMPFLAG) $(INCL) sigma_main.cpp

$(SIGMA_BOOTSTRAP_MPI_OBJS): $(SIGMA_BOOTSTRAP_MPI_SRC)
	$(MPICC) $(CFLAGS) $(SIGMA_BOOTSTRAP_MPI_SRC)

$(SIGMA_JACKKNIFE_MPI_OBJS): $(SIGMA_JACKKNIFE_MPI_SRC)
	$(MPICC) $(CFLAGS) $(SIGMA_JACKKNIFE_MPI_SRC)

