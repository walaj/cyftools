#CC = /opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12
CC=g++

## Required libraries
CEREAL = -I$(HOME)/software/cereal-1.3.2/include ## header only

## Required header-only libraries for UMAP functionality (all are required)
## Could theoretically have KNN only if want KNN graphs but not umap
UMAP =  -I$(HOME)/git/umappp/include -I$(HOME)/git/CppIrlba/include -DHAVE_UMAPPP -DKNNCOLLE_NO_KMKNN -DKNNCOLLE_NO_HNSW -DKNNCOLLE_NO_ANNOY
EIGEN = -I$(HOME)/software/eigen-3.4.0 -DHAVE_EIGEN
KNN = -I$(HOME)/git/knncolle/include -DHAVE_KNNCOLLE

## Very optional header-only libraries for alternative neighor finding algorithms (e.g. in UMAP, etc)
#HNSW = -I$(HOME)/git/hnswlib
KMEANS = -I$(HOME)/git/CppKmeans/include -DHAVE_KMEANS
#ANNOY = -I$(HOME)/git/annoy

## Optional header/lib library for spatial latent dirichlet allocation
LDA = -I$(HOME)/git/supervised-lda/include -DHAVE_LDAPLUSPLUS
LDALIB = -L $(HOME)/git/supervised-lda/lib -lldaplusplus -Wl,-rpath,$(HOME)/git/supervised-lda/lib

## Optional library for HDF5 output
HD5 = -I$(HOME)/software/hdf5-1.14.1/include -DHAVE_HDF5
HD5LIB = -L $(HOME)/software/hdf5-1.14.1/lib -lhdf5_cpp -lhdf5

## Optional header/liblibrary for image-based functions
TIFF = -I$(HOME)/software/libtiff-git/include -DHAVE_TIFFLIB
JPEG = -L/opt/homebrew/opt/jpeg/lib -L/opt/homebrew/Cellar/libpng/1.6.39/lib
TIFFLD = ${JPEG} -llzma $(HOME)/software/libtiff-git/lib/libtiff.a -ljpeg -lz -lpng

## Optional library for KD-tree functions (radial proximity)
KDTREE = -I$(HOME)/git/KDTree -DHAVE_KDTREE
KDLIB = -L $(HOME)/git/KDTree/lib -lKDTree -Wl,-rpath,$(HOME)/git/KDTree/lib

# cysift specific parameter
#USE64BIT=-DUSE_64_BIT

## set the OMP location
ODIR := /opt/homebrew/opt/libomp/include
ifeq ($(wildcard $(ODIR)),)  # If dir does not exist, probably on HMS server
    OMP := -I$(HOME)/software/openmpi-4.1.5/include
    OMPL := -L$(HOME)/software/openmpi-4.1.5/lib -lmpi

	#CAIRO -- module load cairo
    CAIRO=-DHAVE_CAIRO
    #CGAL=-I$(HOME)/software/CGAL-5.5.2/include -DHAVE_CGAL
    CAIROLIB = -lcairo
    CGALLIB = -lgmp
    #BOOST - module load boost

else

    ## Optional libaries for OMP multithreading support
    OMP := -I/opt/homebrew/opt/libomp/include
    OMPL = -L/opt/homebrew/opt/libomp/lib -lomp

    ## Optional library for PDF support
    CAIRO = -I/opt/homebrew/opt/cairo/include -DHAVE_CAIRO
    CAIROLIB = -L /opt/homebrew/opt/cairo/lib -lcairo

    ## Optional header-only library for Voronoi diagram support
    CGAL = -I/opt/homebrew/opt/cgal/include -DHAVE_CGAL

    ## Optional library for better hashing
    BOOST = -I/opt/homebrew/opt/boost/include -DHAVE_BOOST

    #ML includes
    #MLPACK = -I/opt/homebrew/Cellar/mlpack/4.2.0/include -DHAVE_MLPACK
    #ARMADILLO = -I/opt/homebrew/Cellar/armadillo/12.4.1/include
    #ENSMALLEN = -I/opt/homebrew/Cellar/ensmallen/2.19.1/include

    #ML libs
    ###MLPACKL = -L/opt/homebrew/Cellar/mlpack/4.2.0/lib -Wl,-rpath
    #ARMADILLOL = -L/opt/homebrew/Cellar/armadillo/12.4.1/lib -larmadillo -Wl,-rpath,/opt/homebrew/Cellar/armadillo/12.4.1/lib

endif

## check compiler type. If clang then add -Xpreprocessor flag
COMPILER_TYPE := $(shell $(CC) --version)
ifneq (,$(findstring clang,$(COMPILER_TYPE)))
    OPENMP = -Xpreprocessor -fopenmp
else
    OPENMP = -fopenmp
endif

CFLAGS = -O2 -g -std=c++17 $(USE64BIT) $(OPENMP) -I. -I.. \
          $(TIFF) $(KNN) $(UMAP) $(IRLBA) $(EIGEN) $(KMEANS) \
          $(OMP) $(CEREAL) $(LDA) $(HD5) $(KDTREE) $(ENSMALLEN) \
	  $(ARMADILLO) $(MLPACK) $(CAIRO) $(CGAL) $(BOOST)

LDFLAGS = $(OMPL) $(LDALIB) $(HD5LIB) $(KDLIB) ${TIFFLD} $(ARMADILLOL) \
          $(CAIROLIB) $(CGALLIB)

# Specify the source files
SRCS = cysift.cpp \
       cell_table.cpp \
       cell_table_graph.cpp \
       cell_table_lda.cpp \
       polygon.cpp \
       cell_header.cpp \
       cell_graph.cpp \
       cell_flag.cpp \
       cell_utils.cpp \
       cell_processor.cpp \
       cell_row.cpp \
       tiff_reader.cpp \
       tiff_writer.cpp \
       tiff_header.cpp \
       tiff_utils.cpp \
       tiff_ifd.cpp \
       tiff_image.cpp \
       tiff_cp.cpp \
       color_map.cpp

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
