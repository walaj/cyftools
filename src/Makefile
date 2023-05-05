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

## clang
OPENMP = -Xpreprocessor -fopenmp
OMP = -I/opt/homebrew/opt/libomp/include
OMPL = -L/opt/homebrew/opt/libomp/lib -lomp

## GCC
#OPENMP = -fopenmp

BASE=/Users/jeremiahwala/opt/miniconda3

CFLAGS = -O2 -g -std=c++17 $(OPENMP) -I. -I.. $(NETCDF) $(KNN) $(UMAP) $(ARAND) $(IRLBA) $(HNSW) $(EIGEN) $(KMEANS) $(ANNOY) $(OMP)
#CFLAGS = -g -std=c++17 -I. -I.. $(NETCDF) $(KNN) $(UMAP) $(ARAND) $(IRLBA) $(HNSW) $(EIGEN) $(KMEANS) $(ANNOY)

LDFLAGS = $(OMPL)

# Specify the source files
SRCS = cysift.cpp cell_table.cpp polygon.cpp cell_header.cpp cell_graph.cpp cell_flag.cpp cell_utils.cpp cell_processor.cpp

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
