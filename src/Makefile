#CC = /opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12
CC=g++

KNN = -I$(HOME)/git/knncolle/include
UMAP =  -I$(HOME)/git/umappp/include
ARAND = -I/$(HOME)/git/aarand/include
IRLBA = -I$(HOME)/git/CppIrlba/include
HNSW = -I$(HOME)/git/hnswlib
EIGEN = -I$(HOME)/software/eigen-3.4.0
KMEANS = -I$(HOME)/git/CppKmeans/include
ANNOY = -I$(HOME)/git/annoy
CEREAL = -I$(HOME)/software/cereal-1.3.2/include

## set the OMP location
ODIR := /opt/homebrew/opt/libomp/include
ifeq ($(wildcard $(ODIR)),)  # If directory does not exist
    OMP := -I$(HOME)/software/openmpi-4.1.5/include
    OMPL := -L$(HOME)/software/openmpi-4.1.5/lib -lmpi
else
    OMP := -I/opt/homebrew/opt/libomp/include
    OMPL = -L/opt/homebrew/opt/libomp/lib -lomp
endif

## check compiler type. If clang then add -Xpreprocessor flag
COMPILER_TYPE := $(shell $(CC) --version)
ifneq (,$(findstring clang,$(COMPILER_TYPE)))
    OPENMP = -Xpreprocessor -fopenmp
else
    OPENMP = -fopenmp
endif

CFLAGS = -O2 -g -std=c++17 $(OPENMP) -I. -I.. $(NETCDF) $(KNN) $(UMAP) $(ARAND) $(IRLBA) $(EIGEN) $(KMEANS) $(OMP) $(CEREAL) -DKNNCOLLE_NO_HNSW -DKNNCOLLE_NO_ANNOY
LDFLAGS = $(OMPL)

# Specify the source files
SRCS = cysift.cpp cell_table.cpp polygon.cpp cell_header2.cpp cell_graph.cpp cell_flag.cpp cell_utils.cpp cell_processor.cpp cell_row.cpp

# Specify the object files
OBJS = $(SRCS:.cpp=.o)

# Specify the executable file
EXEC = cysift

# The default target
all: $(EXEC)

# Build the executable from the object files
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile the source files into object files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

# Clean up
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean
