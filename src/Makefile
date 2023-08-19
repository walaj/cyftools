#CC = /opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12
CC=g++

MYHOME=$(HOME)
BREWHOME=/opt/homebrew/opt
#USE64BIT=-DUSE_64_BIT

######################################
####### GITHUB SUBCLONE REPOS ########
######################################

## Required libraries
CEREAL = -I../external/cereal/include ## header only

## Required header-only libraries for UMAP functionality (all are required)
## Could theoretically have KNN only if want KNN graphs but not umap
UMAP =  -I../external/umappp/include -I../external/CppIrlba/include -DHAVE_UMAPPP -DKNNCOLLE_NO_KMKNN -DKNNCOLLE_NO_HNSW -DKNNCOLLE_NO_ANNOY
KNN = -I../external/knncolle/include -DHAVE_KNNCOLLE
EIGEN = -I../external/eigen -DHAVE_EIGEN -Wno-deprecated-declarations

## Very optional header-only libraries for alternative neighor finding algorithms (e.g. in UMAP, etc)
KMEANS = -I../external/CppKmeans/include -DHAVE_KMEANS

## Optional library for KD-tree functions (radial proximity)
KDTREE = -I../external/KDTree -DHAVE_KDTREE

######################################
####### OTHER OPTIONAL LIBARARIES ####
######################################

## Optional header/lib library for spatial latent dirichlet allocation
LDA_INCDIR=$(MYHOME)/git/supervised-lda/included
ifneq ($(wildcard $(LDA_INCDIR)/*),)
    LDA=-I$(LDA_INCDIR) -DHAVE_LDAPLUSPLUS
    LDALIB=-L$(MYHOME)/git/supervised-lda/lib -lldaplusplus -Wl,-rpath,$(MYHOME)/git/supervised-lda/lib
else
    $(warning ******************************************************************)
    $(warning Warning: Directory $(LDA_INCDIR) does not exist or is empty. No LDA functionality will be included)
    $(warning To fix, please compile and include/link https://github.com/angeloskath/supervised-lda)
    $(warning ******************************************************************)
endif


####### HMS SERVER SPECIFIC ######
ifneq ($(wildcard $(BREWDIR)),)  # If dir does not exist, probably on HMS server
    OMP = -I$(MYHOME)/software/openmpi-4.1.5/include
    OMPL= -L$(MYHOME)/software/openmpi-4.1.5/lib -lmpi

	#CAIRO -- module load cairo
    CAIRO=-DHAVE_CAIRO
    #CGAL=-I$(MYHOME)/software/CGAL-5.5.2/include -DHAVE_CGAL
    CAIROLIB = -lcairo
    CGALLIB = -lgmp
    #BOOST - module load boost

    ## Optional library for HDF5 output
    # HD5_INCDIR=$(MYHOME)/software/hdf5-1.14.1/included
    # ifneq ($(wildcard $(HD5_INCDIR)/*),)
    #     HD5=-I$(HD5_INCDIR) -DHAVE_HDF5
    #     HD5LIB = -L $(MYHOME)/software/hdf5-1.14.1/lib -lhdf5_cpp -lhdf5
    # else
    #     $(warning ******************************************************************)
    #     $(warning Warning: Directory $(HD5_INCDIR) does not exist or is empty. Output to HDF5 will not be included)
    #     $(warning To fix, please compile and include/link to HDF5 library)
    #     $(warning ******************************************************************)
    # endif    

    ## Optional header/liblibrary for image-based functions
    # TIFF_INCDIR=$(MYHOME)/software/libtiff-git/includedsd
    # ifneq ($(wildcard $(TIFF_INCDIR)/*),) 
    #     TIFF = -I$(TIFF_INCDIR) -DHAVE_TIFFLIB
    #     JPEG = -L/opt/homebrew/opt/jpeg/lib -L/opt/homebrew/Cellar/libpng/1.6.39/lib	
    #     TIFFLD = ${JPEG} -llzma $(MYHOME)/software/libtiff-git/lib/libtiff.a -ljpeg -lz -lpng
    # else
    #     $(warning ******************************************************************)
    #     $(warning Warning: Directory $(TIFF_INCDIR) does not exist or is empty. libtiff will not be included)
    #     $(warning To fix, please compile and include/link to libtif and jpeg libraries)
    #     $(warning ******************************************************************)
    # endif

####### MAC SPECIFIC #######
else

    ## Optional header/liblibrary for image-based functions
    TIFF_INCDIR=$(BREWHOME)/libtiff/include
    ifneq ($(wildcard $(TIFF_INCDIR)/*),) 
        TIFF = -I$(TIFF_INCDIR) -DHAVE_TIFFLIB
        JPEG_INCDIR=$(BREWHOME)/libjpeg
        ifneq ($(wildcard $(JPEG_INCDIR)/*),)
           JPEG = -L$(BREWHOME)/libjpeg/lib -L$(BREWHOME)/libpng/lib -ljpeg -lpng
        else
            $(warning ******************************************************************)
            $(warning Warning: Directory $(JPEG_INCDIR) does not exist or is empty. jpeg support will not be included)
            $(warning To fix, please compile and include/link to jpeg library)
            $(warning ******************************************************************)
        endif
        #TIFFLD = ${JPEG} -llzma $(MYHOME)/software/libtiff-git/lib/libtiff.a -ljpeg -lz -lpng
        TIFFLD = ${JPEG} -llzma -lz -L$(BREWHOME)/libtiff/lib -ltiff
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(TIFF_INCDIR) does not exist or is empty. libtiff will not be included)
        $(warning To fix, please compile and include/link to libtif and jpeg libraries)
        $(warning ******************************************************************)
    endif

    ## Optional library for HDF5 output
    HD5_INCDIR=$(BREWHOME)/hdf5/include
    ifneq ($(wildcard $(HD5_INCDIR)/*),)
        HD5=-I$(HD5_INCDIR) -DHAVE_HDF5
        HD5LIB = -L $(BREWHOME)/hdf5/lib -lhdf5_cpp -lhdf5
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(HD5_INCDIR) does not exist or is empty. Output to HDF5 will not be included)
        $(warning To fix, please compile and include/link to HDF5 library)
        $(warning ******************************************************************)
    endif

    ## Optional libaries for OMP multithreading support
    OMP_INCDIR=$(BREWHOME)/libomp/include
    ifneq ($(wildcard $(OMP_INCDIR)/*),) 
        OMP = -I$(OMP_INCDIR)
        OMPL= -L$(BREWHOME)/libomp/lib -lomp
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(OMP_INCDIR) does not exist or is empty. OMP will not be included)
        $(warning To fix, please include OMP libraries for multithreading support)
        $(warning ******************************************************************)
    endif

    ## Optional library for PDF/PNG support
    CAIRO_INCDIR=$(BREWHOME)/cairo/include
    ifneq ($(wildcard $(CAIRO_INCDIR)/*),) 
        CAIRO = -I$(CAIRO_INCDIR) -DHAVE_CAIROD
        CAIROLIB = -L$(BREWHOME)/cairo/lib -lcairo
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(CAIRO_INCDIR) does not exist or is empty. Cairo will not be included)
        $(warning To fix, please include/link Cairo libraries)
        $(warning ******************************************************************)
    endif

    ## Optional header library for Voronoi support
    CGAL_INCDIR=$(BREWHOME)/cgal/include
    ifneq ($(wildcard $(CGAL_INCDIR)/*),) 
        CGAL = -I$(CGAL_INCDIR) -DHAVE_CGAL
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(CGAL_INCDIR) does not exist or is empty. CGAL will not be included)
        $(warning To fix, please include CGAL libraries)
        $(warning ******************************************************************)
    endif

    ## Optional library for better hashing
    BOOST_INCDIR=$(BREWHOME)/boost/include
    ifneq ($(wildcard $(BOOST_INCDIR)/*),) 
        BOOST = -I$(BOOST_INCDIR) -DHAVE_BOOST
    else
        $(warning ******************************************************************)
        $(warning Warning: Directory $(BOOST_INCDIR) does not exist or is empty. boost will not be included)
        $(warning To fix, please include Boost libraries)
        $(warning ******************************************************************)
    endif

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
