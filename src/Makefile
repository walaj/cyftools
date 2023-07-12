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
LDA = -I$(HOME)/git/supervised-lda/include -DHAVE_LDAPLUSPLUS
HD5 = -I$(HOME)/software/hdf5-1.14.1/include
TIFF = -I$(HOME)/software/libtiff-git/include
KDTREE = -I$(HOME)/git/KDTree


LDALIB = -L $(HOME)/git/supervised-lda/lib -lldaplusplus -Wl,-rpath,$(HOME)/git/supervised-lda/lib

HD5LIB = -L $(HOME)/software/hdf5-1.14.1/lib -lhdf5_cpp -lhdf5
KDLIB = -L $(HOME)/git/KDTree/lib -lKDTree -Wl,-rpath,$(HOME)/git/KDTree/lib
TIFFLD = ${JPEG} -llzma $(HOME)/software/libtiff-git/lib/libtiff.a -ljpeg -lz -lpng

# cysift specific parameter
#USE64BIT=-DUSE_64_BIT

## set the OMP location
ODIR := /opt/homebrew/opt/libomp/include
ifeq ($(wildcard $(ODIR)),)  # If dir does not exist, probably on HMS server
    OMP := -I$(HOME)/software/openmpi-4.1.5/include
    OMPL := -L$(HOME)/software/openmpi-4.1.5/lib -lmpi

	#CAIRO -- module load cairo
    CGAL=-I$(HOME)/software/CGAL-5.5.2/include
    CAIROLIB = -lcairo
    CGALLIB = -lgmp
    #BOOST - module load boost

else
    OMP := -I/opt/homebrew/opt/libomp/include
    OMPL = -L/opt/homebrew/opt/libomp/lib -lomp

    CAIRO = -I/opt/homebrew/opt/cairo/include
    CAIROLIB = -L /opt/homebrew/opt/cairo/lib -lcairo
    CGAL = -I/opt/homebrew/opt/cgal/include
    BOOST = -I/opt/homebrew/opt/boost/include
    JPEG = -L/opt/homebrew/opt/jpeg/lib -L/opt/homebrew/Cellar/libpng/1.6.39/lib

    #ML includes
    MLPACK = -I/opt/homebrew/Cellar/mlpack/4.2.0/include -DHAVE_MLPACK
    ARMADILLO = -I/opt/homebrew/Cellar/armadillo/12.4.1/include
    ENSMALLEN = -I/opt/homebrew/Cellar/ensmallen/2.19.1/include

    #ML libs
    #MLPACKL = -L/opt/homebrew/Cellar/mlpack/4.2.0/lib -Wl,-rpath
    ARMADILLOL = -L/opt/homebrew/Cellar/armadillo/12.4.1/lib -larmadillo -Wl,-rpath,/opt/homebrew/Cellar/armadillo/12.4.1/lib


endif

## check compiler type. If clang then add -Xpreprocessor flag
COMPILER_TYPE := $(shell $(CC) --version)
ifneq (,$(findstring clang,$(COMPILER_TYPE)))
    OPENMP = -Xpreprocessor -fopenmp
else
    OPENMP = -fopenmp
endif

CFLAGS = -O2 -g -std=c++17 $(USE64BIT) $(OPENMP) -I. -I.. $(TIFF) $(KNN) $(UMAP) $(ARAND) $(IRLBA) $(EIGEN) $(KMEANS) $(OMP) $(CEREAL) -DKNNCOLLE_NO_HNSW -DKNNCOLLE_NO_ANNOY $(LDA) $(HD5) $(KDTREE) $(ENSMALLEN) $(ARMADILLO) $(MLPACK) $(CAIRO) $(CGAL) $(BOOST)
LDFLAGS = $(OMPL) $(LDALIB) $(HD5LIB) $(KDLIB) ${TIFFLD} $(ARMADILLOL) $(CAIROLIB) $(CGALLIB)

# Specify the source files
SRCS = cysift.cpp cell_table.cpp polygon.cpp cell_header.cpp cell_graph.cpp cell_flag.cpp cell_utils.cpp cell_processor.cpp cell_row.cpp tiff_reader.cpp tiff_writer.cpp tiff_header.cpp tiff_utils.cpp tiff_ifd.cpp tiff_image.cpp tiff_cp.cpp

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
